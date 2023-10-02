#!/usr/bin/env bash

redo-ifchange all.rid-rcnm.dct.do {linkoutenz,pages,bairoch}.rid-rcnm.dct
lists2list.py {linkoutenz,pages,bairoch}.rid-rcnm.dct > "$3"
