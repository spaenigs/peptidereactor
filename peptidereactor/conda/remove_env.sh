#!/usr/bin/env bash

read -p "Delete environment 'peptidereactor-env'? (y/n)" -n 1 -r
echo

if [[ $REPLY =~ ^[Yy]$ ]]
then
  conda remove --name peptidereactor-env --all
fi