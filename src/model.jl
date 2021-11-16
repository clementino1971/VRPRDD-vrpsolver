function build_model(data::DataVRPRDD, app)

   A = arcs(data) 
   V = [i for i in 1:n(data)] 
   Q = veh_capacity(data)

   # Formulation
   vrprd = VrpModel()
   @variable(vrprd.formulation, x[a in A], Int)
   @variable(vrprd.formulation, y[a in A], Int)

   @objective(vrprd.formulation, Min, sum(c(data,a) * x[a] for a in A))
   #@objective(vrprd.formulation, Min, (sum(c(data,a) * x[a] for a in A)) + (sum((c(data,a)+penalty(data,a[2])) * y[a]) for a in A))
   @constraint(vrprd.formulation, indeg[i in V], sum((x[a]+y[a]) for a in A if a[2] == i) == 1.0)

   println(vrprd.formulation)

   function buildgraph(release_date::Int)
      v_source = v_sink = 0
     
      V1 = [0]
      V2 = Int64[]
      for i in 1:n(data)
         if(release(data, i) == release_date)
            push!(V1,i)
         elseif (release(data, i) < release_date && (release_date + t(data, (i, 0)) <= due(data,i)))
            push!(V1,i)
         elseif (release(data, i) < release_date && (release_date + t(data, (i, 0)) <= deadline(data,i)))
            push!(V2,i)
         end   
      end

      println(V1)
      println(V2)
      
      L, U = 0, upperBoundNbVehicles(data) # multiplicity

      G = VrpGraph(vrprd, vcat(V1,V2), v_source, v_sink, (L, U))

      if app["enable_cap_res"]
         cap_res_id = add_resource!(G, main = true)
      end
      time_res_id = add_resource!(G, main = true)
      
      for v in vcat(V1,V2)
         if app["enable_cap_res"]
            set_resource_bounds!(G, v, cap_res_id, 0, Q)
         end

        #set_resource_bounds!(G, v, time_res_id, release_date, u(data, v))
      end

      for (i,j) in A
         if(i in V1 && j in V1)
            arc_due_id = add_arc!(G, i, j)
            add_arc_var_mapping!(G, arc_due_id, x[(i,j)])
            set_arc_resource_bounds!(G, arc_due_id, time_res_id, release_date, due(data,j))
            set_arc_consumption!(G, arc_due_id, time_res_id, t(data, (i, j)))
            
            if(j != 0)
               arc_dead_id = add_arc!(G, i, j)
               add_arc_var_mapping!(G, arc_dead_id, y[(i,j)])
               set_arc_resource_bounds!(G, arc_dead_id, time_res_id, due(data,j), deadline(data,j))
               set_arc_consumption!(G, arc_dead_id, time_res_id, t(data, (i, j)))
            end

            if app["enable_cap_res"]
               set_arc_consumption!(G, arc_due_id, cap_res_id, d(data, j))
               set_arc_consumption!(G, arc_dead_id, cap_res_id, d(data, j))
            end

         elseif (j in V2 && i in vcat(V1,V2))
            arc_dead_id = add_arc!(G, i, j)
            add_arc_var_mapping!(G, arc_dead_id, y[(i,j)])
            set_arc_resource_bounds!(G, arc_dead_id, time_res_id, due(data,j), deadline(data,j))


            if app["enable_cap_res"]
               set_arc_consumption!(G, arc_dead_id, cap_res_id, d(data, j))
            end

            set_arc_consumption!(G, arc_dead_id, time_res_id, t(data, (i, j)))
         end
      end

      return G,vcat(V1,V2)
   end
   
   all_release_dates = []
   for i in 1:n(data)
      if(! (release(data, i) in all_release_dates) )
         push!(all_release_dates, release(data, i))
      end
   end

   #println(all_release_dates)

   graphs, packing_set_vertex = [], [[] for v in 1:n(data)]
   k = 1
   for rd in all_release_dates
      G,V1 = buildgraph(Int(rd))
      add_graph!(vrprd, G)   
      println(G)
      push!(graphs,G)
      #println(V1, packing_set_vertex)
      for v in V1
         if(v == 0)
            continue
         end
         push!(packing_set_vertex[v],k)
      end
      k+=1
   end

   #we put here a arc packing sets maybe 
   set_vertex_packing_sets!(vrprd, [[(graphs[k],v) for k in packing_set_vertex[v]]  for v in 1:n(data)])
   add_capacity_cut_separator!(vrprd, [ ( [(graphs[k],v) for k in packing_set_vertex[v]], Float64(d(data, v))) for v in 1:n(data)], Float64(Q))

   set_branching_priority!(vrprd, "x", 1)

   return (vrprd, x, y)
end
