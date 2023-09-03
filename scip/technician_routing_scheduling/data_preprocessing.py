# %% [markdown]
# Code taken from the [Technician Routing and Scheduling model](https://github.com/Gurobi/modeling-examples/blob/master/technician_routing_scheduling/technician_routing_scheduling.ipynb)
# from Gurobi.
# %%
import sys
from collections import defaultdict
import xlrd

class Technician():
    def __init__(self, name, cap, depot):
        self.name = name
        self.cap = cap
        self.depot = depot

    def __str__(self):
        return f"Technician: {self.name}\n  Capacity: {self.cap}\n  Depot: {self.depot}"

class Job():
    def __init__(self, name, priority, duration, coveredBy):
        self.name = name
        self.priority = priority
        self.duration = duration
        self.coveredBy = coveredBy

    def __str__(self):
        about = f"Job: {self.name}\n  Priority: {self.priority}\n  Duration: {self.duration}\n  Covered by: "
        about += ", ".join([t.name for t in self.coveredBy])
        return about

class Customer():
    def __init__(self, name, loc, job, tStart, tEnd, tDue):
        self.name = name
        self.loc = loc
        self.job = job
        self.tStart = tStart
        self.tEnd = tEnd
        self.tDue = tDue

    def __str__(self):
        coveredBy = ", ".join([t.name for t in self.job.coveredBy])
        return f"Customer: {self.name}\n  Location: {self.loc}\n  Job: {self.job.name}\n  Priority: {self.job.priority}\n  Duration: {self.job.duration}\n  Covered by: {coveredBy}\n  Start time: {self.tStart}\n  End time: {self.tEnd}\n  Due time: {self.tDue}"

def read_and_process_data(excel_file='data-Sce0.xls'):
    # Open Excel workbook
    wb = xlrd.open_workbook('data\\'+excel_file)

    # Read technician data
    ws = wb.sheet_by_name('Technicians')
    technicians = []
    for i,t in enumerate(ws.col_values(0)[3:]):
        # Create Technician object
        thisTech = Technician(*ws.row_values(3+i)[:3])
        technicians.append(thisTech)

    # Read job data
    jobs = []
    for j,b in enumerate(ws.row_values(0)[3:]):
        coveredBy = [t for i,t in enumerate(technicians) if ws.cell_value(3+i,3+j) == 1]
        # Create Job object
        thisJob = Job(*ws.col_values(3+j)[:3], coveredBy)
        jobs.append(thisJob)

    # Read location data
    ws = wb.sheet_by_name('Locations')
    locations = ws.col_values(0)[1:]
    dist = {(l, l) : 0 for l in locations}
    for i,l1 in enumerate(locations):
        for j,l2 in enumerate(locations):
            if i < j:
                dist[l1,l2] = ws.cell_value(1+i, 1+j)
                dist[l2,l1] = dist[l1,l2]

    # Read customer data
    ws = wb.sheet_by_name('Customers')
    customers = []
    for i,c in enumerate(ws.col_values(0)[1:]):
        for b in jobs:
            if b.name == ws.cell_value(1+i, 2):
                # Create Customer object using corresponding Job object
                rowVals = ws.row_values(1+i)
                #print(rowVals)
                thisCustomer = Customer(*rowVals[:2], b, *rowVals[3:])
                customers.append(thisCustomer)
                break
    
    return technicians, customers, dist

#%%
if __name__=='__main__':
    technicians, customers, dist = read_and_process_data(excel_file='')