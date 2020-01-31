FROM jupyter/scipy-notebook

EXPOSE 8888

COPY ./src /home/jovyan/work

CMD start-notebook.sh --NotebookApp.password='sha1:97756b342378:9974c5dda1cdb015077383fde55f0f8bd4b44621'
