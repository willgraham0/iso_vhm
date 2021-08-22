FROM python:3.8.1-buster

ENV DISPLAY=${DISPLAY}

RUN apt-get update && \
    apt-get install -y build-essential libgtk-3-dev && \
    pip install --no-cache-dir --disable-pip-version-check -U -f https://extras.wxpython.org/wxPython4/extras/linux/gtk3/debian-9 wxPython

WORKDIR /opt/iso-vhm

COPY . .

RUN adduser --system --uid 1000 --group gui-user && \
    chown -R gui-user:gui-user /opt/iso-vhm

USER gui-user

RUN pip install --no-cache-dir --disable-pip-version-check -r requirements.txt

CMD python main.py
