FROM python:3.7

WORKDIR /usr/src/app

COPY resources/requirements.txt ./
COPY resources/entrypoint.sh ./

RUN pip install -r requirements.txt

VOLUME /wd
WORKDIR /wd

ENV PYTHONPATH /wd

RUN mkdir /.streamlit
RUN chmod 777 /.streamlit

ENTRYPOINT ["sh", "/wd/peptidereactor-vis/resources/entrypoint.sh"]