MODULE h5_wrapper
  USE HDF5
  USE constants
  IMPLICIT NONE

  INTEGER(HID_T), PRIVATE   :: file_id
  INTEGER(HID_T), PRIVATE   :: group_id
  INTEGER(HID_T), PRIVATE   :: dset_id
  INTEGER(HID_T), PRIVATE   :: dspace_id
  INTEGER, PRIVATE          :: rank
  INTEGER, PRIVATE          :: error
  INTEGER(HSIZE_T), DIMENSION(1), PRIVATE   :: data_dim
  INTEGER(HSIZE_T), DIMENSION(1:7), PRIVATE :: data_dims
  INTEGER, PARAMETER, PRIVATE               :: delay = 50

CONTAINS
  SUBROUTINE create_h5(file_name)
    CHARACTER(LEN=*), INTENT(IN)  :: file_name

    ! Initialize HDF5_FORTRAN interface.
    CALL h5open_f(error)
    ! Create a new file using default properties.
    CALL h5fcreate_f(TRIM(file_name), H5F_ACC_TRUNC_F, file_id, error)
    ! Close the file.
    CALL h5fclose_f(file_id, error)
    ! Close HDF5_FORTRAN interface.
    CALL h5close_f(error)

    CALL SLEEPQQ (delay)
  END SUBROUTINE create_h5

  SUBROUTINE write_h5_d1(file_name, dset_name, dset)
    CHARACTER(LEN=*), INTENT(IN)  :: file_name
    CHARACTER(LEN=*), INTENT(IN)  :: dset_name
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)     :: dset

    rank = 1
    data_dim = SHAPE(dset)

    ! Initialize HDF5_FORTRAN interface.
    CALL h5open_f(error)
    ! Open the file.
    CALL h5fopen_f(TRIM(file_name), H5F_ACC_RDWR_F, file_id, error)
    ! Create the dataspace.
    CALL h5screate_simple_f(rank, data_dim, dspace_id, error)
    ! Create the dataset with default properties.
    CALL h5dcreate_f(file_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    ! Write the dataset.
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dset, data_dim, error)
    ! Terminate access to the data space.
    CALL h5sclose_f(dspace_id, error)
    ! End access to the dataset and release resources used by it.
    CALL h5dclose_f(dset_id, error)
    ! Close the file.
    CALL h5fclose_f(file_id, error)
    ! Close HDF5_FORTRAN interface.
    CALL h5close_f(error)

    CALL SLEEPQQ (delay)
  END SUBROUTINE write_h5_d1

  SUBROUTINE write_h5_d3(file_name, dset_name, dset)
    CHARACTER(LEN=*), INTENT(IN)  :: file_name
    CHARACTER(LEN=*), INTENT(IN)  :: dset_name
    REAL(KIND=dp), DIMENSION(:,:,:), INTENT(IN) :: dset

    rank = 3
    data_dims(1:rank) = SHAPE(dset)

    ! Initialize HDF5_FORTRAN interface.
    CALL h5open_f(error)
    ! Open the file.
    CALL h5fopen_f(TRIM(file_name), H5F_ACC_RDWR_F, file_id, error)
    ! Create the dataspace.
    CALL h5screate_simple_f(rank, data_dims(1:rank), dspace_id, error)
    ! Create the dataset with default properties.
    CALL h5dcreate_f(file_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    ! Write the dataset.
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dset, data_dims(1:rank), error)
    ! Terminate access to the data space.
    CALL h5sclose_f(dspace_id, error)
    ! End access to the dataset and release resources used by it.
    CALL h5dclose_f(dset_id, error)
    ! Close the file.
    CALL h5fclose_f(file_id, error)
    ! Close HDF5_FORTRAN interface.
    CALL h5close_f(error)

    CALL SLEEPQQ (delay)
  END SUBROUTINE write_h5_d3

  SUBROUTINE write_h5_d5(file_name, dset_name, dset)
    CHARACTER(LEN=*), INTENT(IN)  :: file_name
    CHARACTER(LEN=*), INTENT(IN)  :: dset_name
    REAL(KIND=dp), DIMENSION(:,:,:,:,:), INTENT(IN)  :: dset

    rank = 5
    data_dims(1:rank) = SHAPE(dset)

    ! Initialize HDF5_FORTRAN interface.
    CALL h5open_f(error)
    ! Open the file.
    CALL h5fopen_f(TRIM(file_name), H5F_ACC_RDWR_F, file_id, error)
    ! Create the dataspace.
    CALL h5screate_simple_f(rank, data_dims(1:rank), dspace_id, error)
    ! Create the dataset with default properties.
    CALL h5dcreate_f(file_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    ! Write the dataset.
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dset, data_dims(1:rank), error)
    ! Terminate access to the data space.
    CALL h5sclose_f(dspace_id, error)
    ! End access to the dataset and release resources used by it.
    CALL h5dclose_f(dset_id, error)
    ! Close the file.
    CALL h5fclose_f(file_id, error)
    ! Close HDF5_FORTRAN interface.
    CALL h5close_f(error)

    CALL SLEEPQQ (delay)
  END SUBROUTINE write_h5_d5

  SUBROUTINE write_h5_z5(file_name, dset_name, dset)
    CHARACTER(LEN=*), INTENT(IN)  :: file_name
    CHARACTER(LEN=*), INTENT(IN)  :: dset_name
    COMPLEX(KIND=dp), DIMENSION(:,:,:,:,:), INTENT(IN)  :: dset

    CALL write_h5_d5(file_name,dset_name // '-re', REAL(dset))
    CALL write_h5_d5(file_name,dset_name // '-im', AIMAG(dset))

  END SUBROUTINE write_h5_z5

END MODULE h5_wrapper
