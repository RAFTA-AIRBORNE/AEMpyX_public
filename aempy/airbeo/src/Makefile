# makefile created by mkmf $id: mkmf,v 18.0 2010/03/02 23:26:08 fms exp $ 
f90=gfortran
f90flags= -O -g 
#-fbounds-check -Wall
ld=gfortran
# ldflags=""


.default:
	-echo $@ does not exist.
exec=airbeo.x
all: exec
cnvrt2_depth.o: ./cnvrt2_depth.f90
	$(f90) $(f90flags) $(otherflags) -c	./cnvrt2_depth.f90
cnvrt2_mpar.o: ./cnvrt2_mpar.f90
	$(f90) $(f90flags) $(otherflags) -c	./cnvrt2_mpar.f90
cnvrt2_xpar.o: ./cnvrt2_xpar.f90
	$(f90) $(f90flags) $(otherflags) -c	./cnvrt2_xpar.f90
cnvrt_bounds.o: ./cnvrt_bounds.f90
	$(f90) $(f90flags) $(otherflags) -c	./cnvrt_bounds.f90
config_id.o: ./config_id.f90
	$(f90) $(f90flags) $(otherflags) -c	./config_id.f90
costrn.o: ./costrn.f90 filt_coef_q.o
	$(f90) $(f90flags) $(otherflags) -c	./costrn.f90
cubder.o: ./cubder.f90
	$(f90) $(f90flags) $(otherflags) -c	./cubder.f90
cubint.o: ./cubint.f90
	$(f90) $(f90flags) $(otherflags) -c	./cubint.f90
cubspl.o: ./cubspl.f90
	$(f90) $(f90flags) $(otherflags) -c	./cubspl.f90
cubval.o: ./cubval.f90
	$(f90) $(f90flags) $(otherflags) -c	./cubval.f90
dcprm_fd.o: ./dcprm_fd.f90
	$(f90) $(f90flags) $(otherflags) -c	./dcprm_fd.f90
dcprm_td.o: ./dcprm_td.f90
	$(f90) $(f90flags) $(otherflags) -c	./dcprm_td.f90
dprod1.o: ./dprod1.f90
	$(f90) $(f90flags) $(otherflags) -c	./dprod1.f90
esvd.o: ./esvd.f90
	$(f90) $(f90flags) $(otherflags) -c	./esvd.f90
filter_coefficients.o: ./filter_coefficients.f90
	$(f90) $(f90flags) $(otherflags) -c	./filter_coefficients.f90
filt_coef_q.o: ./filt_coef_q.f90
	$(f90) $(f90flags) $(otherflags) -c	./filt_coef_q.f90
fold_and_convolve.o: ./fold_and_convolve.f90
	$(f90) $(f90flags) $(otherflags) -c	./fold_and_convolve.f90
forjac.o: ./forjac.f90
	$(f90) $(f90flags) $(otherflags) -c	./forjac.f90
frequency_select.o: ./frequency_select.f90
	$(f90) $(f90flags) $(otherflags) -c	./frequency_select.f90
hsboss_fd.o: ./hsboss_fd.f90
	$(f90) $(f90flags) $(otherflags) -c	./hsboss_fd.f90
hsboss_td.o: ./hsboss_td.f90 frequency_select.o
	$(f90) $(f90flags) $(otherflags) -c	./hsboss_td.f90
hsmd_fd.o: ./hsmd_fd.f90
	$(f90) $(f90flags) $(otherflags) -c	./hsmd_fd.f90
hsmd_hnk.o: ./hsmd_hnk.f90 filt_coef_q.o
	$(f90) $(f90flags) $(otherflags) -c	./hsmd_hnk.f90
hsmd_ker.o: ./hsmd_ker.f90
	$(f90) $(f90flags) $(otherflags) -c	./hsmd_ker.f90
input_data_for_airbeo.o: ./input_data_for_airbeo.f90
	$(f90) $(f90flags) $(otherflags) -c	./input_data_for_airbeo.f90
interv.o: ./interv.f90
	$(f90) $(f90flags) $(otherflags) -c	./interv.f90
linval.o: ./linval.f90
	$(f90) $(f90flags) $(otherflags) -c	./linval.f90
airbeo.o: ./airbeo.f90 input_data_for_airbeo.o
	$(f90) $(f90flags) $(otherflags) -c	./airbeo.f90
nlsq2.o: ./nlsq2.f90
	$(f90) $(f90flags) $(otherflags) -c	./nlsq2.f90
noise_2_sigr.o: ./noise_2_sigr.f90
	$(f90) $(f90flags) $(otherflags) -c	./noise_2_sigr.f90
parameter_sensitivity.o: ./parameter_sensitivity.f90
	$(f90) $(f90flags) $(otherflags) -c	./parameter_sensitivity.f90
set_norm_td.o: ./set_norm_td.f90
	$(f90) $(f90flags) $(otherflags) -c	./set_norm_td.f90
set_source.o: ./set_source.f90
	$(f90) $(f90flags) $(otherflags) -c	./set_source.f90
solve2.o: ./solve2.f90
	$(f90) $(f90flags) $(otherflags) -c	./solve2.f90
tqstrip.o: ./tqstrip.f90
	$(f90) $(f90flags) $(otherflags) -c	./tqstrip.f90
txcmrg.o: ./txcmrg.f90
	$(f90) $(f90flags) $(otherflags) -c	./txcmrg.f90
txcnvd.o: ./txcnvd.f90
	$(f90) $(f90flags) $(otherflags) -c	./txcnvd.f90
txcnvl.o: ./txcnvl.f90
	$(f90) $(f90flags) $(otherflags) -c	./txcnvl.f90
write_fd.o: ./write_fd.f90
	$(f90) $(f90flags) $(otherflags) -c	./write_fd.f90
write_log_file.o: ./write_log_file.f90
	$(f90) $(f90flags) $(otherflags) -c	./write_log_file.f90
write_misfit.o: ./write_misfit.f90
	$(f90) $(f90flags) $(otherflags) -c	./write_misfit.f90
write_model.o: ./write_model.f90
	$(f90) $(f90flags) $(otherflags) -c	./write_model.f90
write_td.o: ./write_td.f90
	$(f90) $(f90flags) $(otherflags) -c	./write_td.f90
wrslv.o: ./wrslv.f90
	$(f90) $(f90flags) $(otherflags) -c	./wrslv.f90
wrslvp.o: ./wrslvp.f90
	$(f90) $(f90flags) $(otherflags) -c	./wrslvp.f90
wrslvs_fd.o: ./wrslvs_fd.f90
	$(f90) $(f90flags) $(otherflags) -c	./wrslvs_fd.f90
wrslv_fd.o: ./wrslv_fd.f90
	$(f90) $(f90flags) $(otherflags) -c	./wrslv_fd.f90
src = ./cnvrt2_xpar.f90 ./hsboss_td.f90 ./cnvrt2_mpar.f90 ./noise_2_sigr.f90 ./write_fd.f90 ./hsmd_hnk.f90 ./cubint.f90 ./tqstrip.f90 ./parameter_sensitivity.f90 ./txcmrg.f90 ./esvd.f90 ./wrslvs_fd.f90 ./write_misfit.f90 ./linval.f90 ./set_norm_td.f90 ./cubspl.f90 ./cnvrt_bounds.f90 ./fold_and_convolve.f90 ./wrslv_fd.f90 ./interv.f90 ./cubval.f90 ./costrn.f90 ./hsmd_ker.f90 ./write_model.f90 ./dcprm_td.f90 ./filter_coefficients.f90 ./input_data_for_airbeo.f90 ./dcprm_fd.f90 ./write_td.f90 ./txcnvd.f90 ./set_source.f90 ./forjac.f90 ./solve2.f90 ./wrslvp.f90 ./dprod1.f90 ./txcnvl.f90 ./airbeo.f90 ./hsboss_fd.f90 ./wrslv.f90 ./write_log_file.f90 ./nlsq2.f90 ./cubder.f90 ./cnvrt2_depth.f90 ./filt_coef_q.f90 ./frequency_select.f90 ./hsmd_fd.f90 ./config_id.f90
obj = cnvrt2_xpar.o hsboss_td.o cnvrt2_mpar.o noise_2_sigr.o write_fd.o hsmd_hnk.o cubint.o tqstrip.o parameter_sensitivity.o txcmrg.o esvd.o wrslvs_fd.o write_misfit.o linval.o set_norm_td.o cubspl.o cnvrt_bounds.o fold_and_convolve.o wrslv_fd.o interv.o cubval.o costrn.o hsmd_ker.o write_model.o dcprm_td.o filter_coefficients.o input_data_for_airbeo.o dcprm_fd.o write_td.o txcnvd.o set_source.o forjac.o solve2.o wrslvp.o dprod1.o txcnvl.o airbeo.o hsboss_fd.o wrslv.o write_log_file.o nlsq2.o cubder.o cnvrt2_depth.o filt_coef_q.o frequency_select.o hsmd_fd.o config_id.o

# exmplist=
# tests:  $(exec) 
#	@for file in $(exmplist); do  cp file;$(exec)  ; done

clean: 
	-rm -f $(tmpfiles) $(obj)   $(exec)  *.mod
neat:
	-rm -f $(tmpfiles)

exec: $(obj) 
	$(ld) $(obj)  -o  $(exec)  $(ldflags)
