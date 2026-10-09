"""Process-wide lock for netCDF reads: netCDF4/HDF5 is not thread-safe, and the server and the
experiment tools run jobs in threads (concurrent opens give "Not a valid ID" / "HDF error")."""
import threading

NC_LOCK = threading.RLock()
