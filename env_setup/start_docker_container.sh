sudo docker run --rm --shm-size=256m -p 8888:8888 -e PASSWORD=huhu \
-e JUPYTER_ENABLE_LAB=yes -e GRANT_SUDO=yes --user root \
-v ${working_dir}/3D_gene_expression_prediction/:/home/jovyan/storage/ \
--name jupyter_lab_novosparc jupyter_lab_novosparc &

# sudo docker run --rm -it --shm-size=256m -p 8888:8888 -e PASSWORD=huhu \
# -e JUPYTER_ENABLE_LAB=yes -e GRANT_SUDO=yes --user root \
# -v ${working_dir}/3D_gene_expression_prediction/:/home/jovyan/storage/ \
# --name jupyter_lab_novosparc jupyter_lab_novosparc &
