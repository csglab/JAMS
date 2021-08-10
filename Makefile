MAKE = make		#change this line if you are using a different GNU make software


dirAffiMx = ./src/affimx
hg19 = ./data/hg19
mmaps = ./data/methyl_maps/hek293

all: CC_AffiMx Extract 

Extract:
	gunzip $(hg19)/*.gz $(mmaps)/*.gz

CC_AffiMx: $(dirAffiMx)/Makefile
	$(MAKE) -C $(dirAffiMx)

clean:
	echo "Done"
