#!/bin/bash

GRP="$1"
newgrp $GRP
export HOME=/home/$GRP/zhoux379
cd
source .bashrc

