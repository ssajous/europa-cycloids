FROM jupyter/scipy-notebook

EXPOSE 8888

COPY ./src /home/jovyan/work

CMD start-notebook.sh