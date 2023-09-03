# %% [markdown]
# One to one translation from the 
# [Technician Routing and Scheduling model](https://github.com/Gurobi/modeling-examples/blob/master/technician_routing_scheduling/technician_routing_scheduling.ipynb) 
# provided by the Gurobi documentation. 
# 
# The following mathematical model is also copied from the Gurobi documentation for describing the code at hand.
#%%
from data_preprocessing import read_and_process_data
from pyscipopt import Model, quicksum

# %%
def trs_model(technicians, customers, dist, time_limit, gap_limit):
    # Taking same data structures as Gurobi
    # Build useful data structures
    K = [k.name for k in technicians]
    C = [j.name for j in customers]
    J = [j.loc for j in customers]
    L = list(set([l[0] for l in dist.keys()]))
    D = list(set([t.depot for t in technicians]))
    cap = {k.name : k.cap for k in technicians}
    loc = {j.name : j.loc for j in customers}
    depot = {k.name : k.depot for k in technicians}
    canCover = {j.name : [k.name for k in j.job.coveredBy] for j in customers}
    dur = {j.name : j.job.duration for j in customers}
    tStart = {j.name : j.tStart for j in customers}
    tEnd = {j.name : j.tEnd for j in customers}
    tDue = {j.name : j.tDue for j in customers}
    priority = {j.name : j.job.priority for j in customers}

    model = Model('Model')
    if time_limit != None:
        model.setParam('limits/time', time_limit)
    if gap_limit != None:
        model.setParam('limits/gap', gap_limit)

    ### Variables ###
    # Customer-technician assignment
    x = {}
    for c in C:
        for k in K:
            x[c,k] = model.addVar(vtype="B", name="x(%s,%s)"%(c,k))    
    # Technician assignment
    u = {}
    for k in K:
        u[k] = model.addVar(vtype="B", name="x(%s)"%(k))
    
    # Edge-route assignment to technician
    y = {}
    for i in L:
        for j in L:
            for k in K:
                y[i,j,k] = model.addVar(vtype="B", name="y(%s,%s,%s)"%(i,j,k))
   
    # Technician cannot leave or return to a depot that is not its base
    for k in technicians:
        for d in D:
            if k.depot != d:
                for i in L:
                    y[i,d,k.name] = model.addVar(ub=0, vtype="B", name="y(%s,%s,%s)"%(i,d,k.name))
                    y[d,i,k.name] = model.addVar(ub=0, vtype="B", name="y(%s,%s,%s)"%(d,i,k.name))
    
    # Start time of service
    t = {}
    for i in L:
        t[i] = model.addVar(ub=600, vtype="I", name="t(%s)"%(i))
    
    # Lateness of service
    z = {}
    for c in C:
        z[c] = model.addVar(vtype="C", name="z(%s)"%(c))
    
    # Artificial variables to correct time window upper and lower limits
    xa, xb = {}, {}
    g = {}
    for j in C:
        xa[j] = model.addVar(vtype="C", name="xa(%s)"%(j))
        xb[j] = model.addVar(vtype="C", name="xb(%s)"%(j))
    
        # Unfilled jobs
        g[j] = model.addVar(vtype="C", name="g(%s)"%(j))
    
    ### Variables ###
    
    ### Constraints ###
    
    for j in C:
        # A technician must be assigned to a job, or a gap is declared (1)
        model.addCons(quicksum(x[j,k] for k in canCover[j]) + g[j] == 1, "assignToJob(%s)"%(j))

        # At most one technician can be assigned to a job (2)
        #model.addCons(quicksum(x[j,k] for k in canCover[j]) <= 1, "assignOne(%s)"%(j))
        model.addCons(quicksum(x[j,k] for k in K) <= 1, "assignOne(%s)"%(j))

    # Technician capacity constraints (3)    
    for k in K:
        model.addCons(quicksum(dur[j]*x[j,k] for j in C) + \
                      quicksum(dist[i,j]*y[i,j,k] for i in L for j in L) <= 
                      (cap[k]*u[k]),
                      "techCapacity(%s)"%(k)
                      )
        
    # Technician tour constraints (4 and 5)
    for k in K:
        for i in C:
            model.addCons(quicksum(y[loc[i],j,k] for j in L) == x[i,k], "techTour1(%s,%s)"%(i,k))
        for j in C:
            model.addCons(quicksum(y[i,loc[j],k] for i in L) == x[j,k], "techTour2(%s,%s)"%(j,k))

    # Same depot constraints (6 and 7)
    for k in K:
        model.addCons((quicksum(y[j,depot[k],k] for j in J) == u[k]), "sameDepot1(%s)"%(k))
        model.addCons((quicksum(y[depot[k],j,k] for j in J) == u[k]), "sameDepot2(%s)"%(k))

    # Temporal constraints (8) for customer locations
    M = {(i,j) : 600 + dur[i] + dist[loc[i], loc[j]] for i in C for j in C}
    for i in C:
        for j in C:
            model.addCons(t[loc[j]] >= t[loc[i]] + dur[i] + dist[loc[i], loc[j]] \
                        - M[i,j]*(1 - quicksum(y[loc[i],loc[j],k] for k in K))                   
                        ,"tempoCustomer(%s,%s)"%(i,j))
    
    # Temporal constraints (8) for depot locations
    M = {(i,j) : 600 + dist[i, loc[j]] for i in D for j in C}        
    for i in D:
        for j in C:
            model.addCons(t[loc[j]] >= t[i] + dist[i, loc[j]] \
                        - M[i,j]*(1 - quicksum(y[i,loc[j],k] for k in K))                   
                        ,"tempoDepot(%s,%s)"%(i,j))

    for j in C:
        # Time window constraints (9 and 10)
        model.addCons(t[loc[j]] + xa[j] >= tStart[j], "timeWinA(%s)"%(j))
        model.addCons(t[loc[j]] - xb[j] <= tEnd[j], "timeWinB(%s)"%(j))

        # Lateness constraint (11)
        model.addCons(z[j] >= t[loc[j]] + dur[j] - tDue[j], "lateness(%s)"%(j))

    ### Constraints ###

    ### Objective Function ###
    M = 6100
    
    model.setObjective(quicksum(priority[j]*z[j] for j in C) + \
                       quicksum( 0.01 * M * priority[j] * (xa[j] + xb[j]) for j in C) + \
                       quicksum( M * priority[j] * g[j] for j in C)
                       ,"minimize")
    #model.setObjective(quicksum(g[j] for j in C), "minimize")
    ### Objective Function ###
    
    model.data = x, u, y, t, z, xa, xb, g
    return model

#%%
if __name__=='__main__':
    technicians, customers, dist = read_and_process_data(excel_file='data-Sce3.xls')
    model = trs_model(technicians, customers, dist, time_limit=60, gap_limit=0.001)
    model_data = model.data
    dvar_x, dvar_u, dvar_y, dvar_t, dvar_z, dvar_xa, dvar_xb, dvar_g = model.data
    model.optimize()

#%%
if model.getStatus=='infeasible':
    print()
else:
    print(model.getObjVal(), model.getPrimalbound(), model.getDualbound(), model.getGap())

#%%
rows = []
for (a,b) in model_data[0]:
    if model.getVal(model_data[0][a,b])>0.5:

# %%
K = [k.name for k in technicians]
C = [j.name for j in customers]
J = [j.loc for j in customers]
L = list(set([l[0] for l in dist.keys()]))
D = list(set([t.depot for t in technicians]))

### Print results
# Assignments
print("")
for j in customers:
    if model.getVal(dvar_g[j.name]) > 0.5:
        jobStr = "Nobody assigned to {} ({}) in {}".format(j.name,j.job.name,j.loc)
    else:
        for k in K:
            if model.getVal(dvar_x[j.name,k]) > 0.5:
                jobStr = "{} assigned to {} ({}) in {}. Start at t={:.2f}.".format(k,j.name,
                                                                                   j.job.name,
                                                                                   j.loc,
                                                                                   model.getVal(dvar_t[j.loc])
                                                                                   )
                if model.getVal(model_data[4][j.name]) > 1e-6:
                    jobStr += " {:.2f} minutes late.".format(model.getVal(dvar_z[j.name]))
                if model.getVal(dvar_xa[j.name]) > 1e-6:
                    jobStr += " Start time corrected by {:.2f} minutes.".format(model.getVal(dvar_xa[j.name]))
                if model.getVal(dvar_xb[j.name]) > 1e-6:
                    jobStr += " End time corrected by {:.2f} minutes.".format(model.getVal(dvar_xb[j.name]))
    print(jobStr)

# Technicians
print("")
for k in technicians:
    if model.getVal(dvar_u[k.name]) > 0.5:
        cur = k.depot
        route = k.depot
        while True:
            for j in customers:
                if model.getVal(dvar_y[cur,j.loc,k.name]) > 0.5:
                    route += " -> {} (dist={}, t={:.2f}, proc={})".format(j.loc, dist[cur,j.loc], model.getVal(dvar_t[j.loc]), j.job.duration)
                    cur = j.loc
            for i in D:
                if model.getVal(dvar_y[cur,i,k.name]) > 0.5:
                    route += " -> {} (dist={})".format(i, dist[cur,i])
                    cur = i
                    break
            if cur == k.depot:
                break
        print("{}'s route: {}".format(k.name, route))
    else:
        print("{} is not used".format(k.name)) 
        
#%%
# Utilization
print("")
for k in K:
    used = capLHS[k].getValue()
    total = cap[k]
    util = used / cap[k] if cap[k] > 0 else 0
    print("{}'s utilization is {:.2%} ({:.2f}/{:.2f})".format(k, util,\
        used, cap[k]))
totUsed = sum(capLHS[k].getValue() for k in K)
totCap = sum(cap[k] for k in K)
totUtil = totUsed / totCap if totCap > 0 else 0
print("Total technician utilization is {:.2%} ({:.2f}/{:.2f})".format(totUtil, totUsed, totCap))