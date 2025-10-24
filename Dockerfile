FROM python:3.10

WORKDIR /usr/src/app

COPY requirements.txt ./
RUN apt-get update
RUN apt-get -y install cmake protobuf-compiler
RUN apt-get install build-essential libgl1-mesa-dev -y
RUN apt install libcairo2-dev libxt-dev libgirepository1.0-dev -y
RUN apt-get install python3-pyqt5 -y
RUN apt install python3-gi python3-gi-cairo gir1.2-gtk-3.0 gir1.2-webkit2-4.1 -y
RUN pip install --no-cache-dir -r requirements.txt
COPY . .

RUN pip install .

EXPOSE 5001
CMD [ "python", "src/kinetic_analysis/app_dash.py"]
