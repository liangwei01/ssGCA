import os
for i in [3,9]:
    i = i*10
    if i != 50:
        with open(f"job_tp_{i}.sh","w") as f:
            f.write(f"#!/bin/bash\n\
#OAR -n sample_tp\n\
#OAR -l /nodes=1/core=8,walltime=20:48:00\n\
#OAR --stdout tp.out\n\
#OAR --stderr tp.err\n\
#OAR --project pr-ai4drug\n\
cd /bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/GCN/katy-gene-networks-main_generating_the_result/leave_one_out\n\
singularity exec --bind /bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/GCN/katy-gene-networks-main_generating_the_result/leave_one_out/:/bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/GCN/katy-gene-networks-main_generating_the_result/leave_one_out /bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/imgs/rbioinfo python /bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/GCN/katy-gene-networks-main_generating_the_result/leave_one_out/get_R2_power_law_and_topology_features.py primaryEver all_sample_k_{i} \n\
singularity exec --bind /bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/GCN/katy-gene-networks-main_generating_the_result/leave_one_out/:/bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/GCN/katy-gene-networks-main_generating_the_result/leave_one_out /bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/imgs/rbioinfo python /bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/GCN/katy-gene-networks-main_generating_the_result/leave_one_out/get_R2_power_law_and_topology_features.py primaryNivo all_sample_k_{i} \n\
singularity exec --bind /bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/GCN/katy-gene-networks-main_generating_the_result/leave_one_out/:/bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/GCN/katy-gene-networks-main_generating_the_result/leave_one_out /bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/imgs/rbioinfo python /bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/GCN/katy-gene-networks-main_generating_the_result/leave_one_out/get_R2_power_law_and_topology_features.py metaEver all_sample_k_{i} \n\
singularity exec --bind /bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/GCN/katy-gene-networks-main_generating_the_result/leave_one_out/:/bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/GCN/katy-gene-networks-main_generating_the_result/leave_one_out /bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/imgs/rbioinfo python /bettik/PROJECTS/pr-ai4drug/yinli/network_analysis/GCN/katy-gene-networks-main_generating_the_result/leave_one_out/get_R2_power_law_and_topology_features.py metaNivo all_sample_k_{i}")
        os.system(f"chmod +x job_tp_{i}.sh")
        #os.system(f"oarsub -S job_tp_{i}.sh")
