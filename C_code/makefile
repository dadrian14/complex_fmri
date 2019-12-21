lib.dll: 
	R CMD SHLIB -o lib.dll mat_vec.c RC_interface.c GLS.c complex_Sig=gen.c \ 
		complex_Sig=sig2I.c complex_smoother.c spatial_smooth.c

clean:
	rm *.o lib.dll

