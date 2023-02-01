include $(shell ap.makesim)

hook = hook/build/libhook.so

par.py:
	touch $@

par.conf: par.py gen
	./gen
conf: par.conf

$(hook): hook/hook.cpp
	make -C hook
conf: $(hook)

gen.copy: gen
	@cp gen gen.copy
conf: gen.copy

clean::
	rm -vf par.make partinit.csv
	rm -f bond_*.vtk
	rm -f gen.copy
