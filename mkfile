<config.mk

results/%.RData:	data/%.RDS
	set -x
	mkdir -p results
	/labs/genut/bin/Rscript ./bin/corto.R \
		--RDS $prereq \
		--transpose TRUE \
		--TFs $TFLIST \
		--nbootstraps $BOOTSTRAPS \
		--pvalue $PVALUE \
		--nthreads $NTHREADS \
		--output $target'.build' \
	&& mv $target'.build' $target

clean:V:
        bin/targets | xargs rm -rf
        rm -rf results/
