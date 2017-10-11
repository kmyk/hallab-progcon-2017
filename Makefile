.PHONY: build run json

build:
	- ln -s $(PWD)/Answer.cpp hpc2017/src/Answer.cpp
	$(MAKE) -C hpc2017 all
run:
	$(MAKE) build
	$(MAKE) -C hpc2017 run
view:
	$(MAKE) build
	$(MAKE) -C hpc2017 json | sed '2!d' > hpc2017/viewer/data.json
	x-www-browser http://localhost:8080/?data=data.json
server:
	# npm install -g http-server
	http-server hpc2017/viewer &
