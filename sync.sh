#!/bin/bash
rsync -avz --delete \
  --exclude 'venv/' \
  --exclude '__pycache__/' \
  --exclude '*.pyc' \
  --exclude '.git/' \
  --exclude 'a.out' \
  --exclude '.gitignore' \
  ./ co24btech11023@nsws.mae.jpan.in:~/canm2