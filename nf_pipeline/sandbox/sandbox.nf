//Default params
params.scoring_processes = 5
params.permute_gene_labels = 'yes'

//Default all metrics off, add to list if user turned them on
def metrics_list = []

params.peripheral = false
if(params.peripheral){
    metrics_list.add('peripheral')
}

params.central = false
if(params.central){
    metrics_list.add('central')
}

params.radial = false
if(params.radial){
    metrics_list.add('radial')
}

params.punctate = false
if(params.punctate){
    metrics_list.add('punctate')
}


metrics_ch = Channel.fromList(metrics_list)

print params.scoring_processes

print params.peripheral
print params.central
print params.radial
print params.punctate
metrics_ch.view()

