.PHONY: build run json

build:
	- ln -s $(PWD)/Answer.cpp hpc2017/src/Answer.cpp
	$(MAKE) -C hpc2017 all
run:
	$(MAKE) build
	$(MAKE) -C hpc2017 run
view:
	$(MAKE) build
	$(MAKE) -C hpc2017 json | sed '2!d;s/.*/localdata = "&";/' > hpc2017/viewer/localdata.js
	x-www-browser file:///$(PWD)/hpc2017/viewer/index.html?localdata=
