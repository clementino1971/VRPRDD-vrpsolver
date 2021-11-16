import Base.show, Base.print

mutable struct Vertex
   id_vertex::Int
   pos_x::Float64
   pos_y::Float64
   demand::Int
   release_date::Float64
   due_date::Float64
   deadline_date::Float64
   penalty::Float64
end

# Directed graph
mutable struct InputGraph
   V′::Array{Vertex} # set of vertices
   A::Array{Tuple{Int64,Int64}} # set of edges
   cost::Dict{Tuple{Int64,Int64},Float64} # cost for each arc
   time::Dict{Tuple{Int64,Int64},Float64} # time for each arc
end

mutable struct DataVRPRDD
   n::Int
   G′::InputGraph
   Q::Float64 # vehicle capacity
   K::Int #num vehicles available
end

# Euclidian distance
function distance(v1::Vertex, v2::Vertex)
   x_sq = (v1.pos_x - v2.pos_x)^2
   y_sq = (v1.pos_y - v2.pos_y)^2
   return floor(sqrt(x_sq + y_sq)*10)/10
end

function readVRPRDDData(path_file::String)
   data = Array{Any,1}()
   dist = Array{Any,1}()
   open(path_file) do file
      for line in eachline(file)
         println(line)
         if findfirst("CUST", line) != nothing || findfirst("VEH", line) != nothing || 
	    findfirst("NUMBER", line) != nothing
            continue
         end
      
         for peaceofdata in split(line)
            push!(data, String(peaceofdata))
         end
      end
   end

   n = div(length(data) - 3, 8) - 1
   K = parse(Int, data[2])
   Q = parse(Float64, data[3])

   println(n, " ", K, " ", Q)

   vertices = Vertex[] 
   for i in 0:n
      
      offset = 3 + i*8
      x = parse(Float64, data[offset + 2])     
      y = parse(Float64, data[offset + 3])     
      d = parse(Int, data[offset + 4])     
      l = parse(Int, data[offset + 5])     
      due = parse(Int, data[offset + 6])
      D = parse(Int, data[offset + 7])
      p = parse(Int, data[offset + 8])     

      push!(vertices, Vertex(i, x, y, d, l, due,D,p))
   end
   
   println(vertices)

   A = Tuple{Int64,Int64}[]
   cost = Dict{Tuple{Int64,Int64},Float64}()
   time = Dict{Tuple{Int64,Int64},Float64}()

   function add_arc!(i, j) 
      push!(A, (i,j)) 
      
      cost[(i,j)] = distance(vertices[i+1], vertices[j+1])  
      time[(i,j)] = cost[(i,j)] 
   end

   for i in 1:n
      #arc from depot
      
      add_arc!(0, i)
      #arc to depot
      add_arc!(i, 0)
      for j in 1:n
         if (i != j) 
            add_arc!(i, j)
         end
      end
   end

   DataVRPRDD(n, InputGraph(vertices, A, cost, time), Q, K)
end

arcs(data::DataVRPRDD) = data.G′.A # return set of arcs
function c(data,a) 
   if !(haskey(data.G′.cost, a)) 
      return Inf
   end
   return data.G′.cost[a] 
end     
function t(data,a) 
   if !(haskey(data.G′.time, a)) 
      return Inf
   end
   return data.G′.time[a] 
end

n(data::DataVRPRDD) = data.n # return number of requests
d(data::DataVRPRDD, i) = data.G′.V′[i+1].demand # return demand of i
release(data::DataVRPRDD, i) = data.G′.V′[i+1].release_date
due(data::DataVRPRDD, i) = data.G′.V′[i+1].due_date
deadline(data::DataVRPRDD, i) = data.G′.V′[i+1].deadline_date
penalty(data::DataVRPRDD, i) = data.G′.V′[i+1].penalty
veh_capacity(data::DataVRPRDD) = Int(data.Q)

function lowerBoundNbVehicles(data::DataVRPRDD) 
   return 1
end

function upperBoundNbVehicles(data::DataVRPRDD) 
   return data.K
end