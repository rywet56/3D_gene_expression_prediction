echo $working_dir

cd "${working_dir}/3D_gene_expression_prediction/docker_files"
sudo docker build ${working_dir}/3D_gene_expression_prediction/docker_files -f \
jupyter_lab_novosparc_dockerfile -t jupyter_lab_novosparc
