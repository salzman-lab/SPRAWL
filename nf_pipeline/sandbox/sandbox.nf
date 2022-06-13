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

metrics_ch.view()

//Default params
params.scoring_processes = 5
params.permute_gene_labels = 'no'
params.shrink_factor = 1

//Output the parameter values to a log file
new FileWriter("outputs/${params.run_name}/parameters.txt").with {
    write(String.format('%tF %<tH:%<tM', java.time.LocalDateTime.now())+"\n")
    for ( e in params ) {
        write("params.${e.key} = ${e.value}\n")
    }
    flush()
} 
    


