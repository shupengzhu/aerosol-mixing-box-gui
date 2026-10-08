# SCRAM box app v2

Box model of the CESM-SCRAM aerosol kernel (latest hetfix code), with a selectable
ISORROPIA version (SCRAM-revised v1 or ISORROPIA II 2.3) and a local web GUI.
No cfg files: the GUI writes a JSON run file that the kernel driver reads.

## Status
- [x] `kernel/build_kernel.sh`: builds the CESM-SCRAM `scram_core` as a standalone
      static library from a snapshot of a CESM-SCRAM source tree; records the kernel
      revision in `build/<variant>/KERNEL_VERSION`; whole-archive link check.
- [ ] JSON box driver (`kernel/driver`)
- [ ] ISORROPIA II 2.3 shim for SCRAM's 11-argument `ISOROPIA2`
- [ ] web GUI

## Build the kernel
    ISO_VARIANT=v1 bash kernel/build_kernel.sh
Host modules replaced by stubs (`kernel/stubs`): `cam_abortutils` (endrun),
`perf_mod` (timers, no-op), `time_manager` (get_nstep, set by the driver).
Compile notes: every object, including Fortran 77, is built with `-fopenmp`
because SCRAM declares ISORROPIA common blocks THREADPRIVATE; link with `-lblas`.
