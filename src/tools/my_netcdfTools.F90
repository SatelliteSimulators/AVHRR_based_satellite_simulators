MODULE my_netcdfTools
  !
  ! Salomon.Eliasson@smhi.se 
  !

  USE netcdf, ONLY : &
       nf90_clobber, &
       nf90_close,   &
       nf90_create,  &
       nf90_def_dim, &
       nf90_def_var, &
       nf90_enddef,  &
       nf90_put_att, &
       nf90_global,  &
       nf90_int,     &
       nf90_noerr,   &
       nf90_put_att, &
       nf90_real,    &
       nf90_strerror,&
       nf90_unlimited


  IMPLICIT NONE

  PUBLIC :: preamble_Netcdf, add_variable_attributes,&
       check          

CONTAINS

  ! ---------------------
  ! NETCDF writing functions
  ! ---------------------

  SUBROUTINE preamble_Netcdf(file_name,  &
       dimtname, dimtid, &
       dim1name, dim1id, dim1len, &
       dim2name, dim2id, dim2len, &
       dim3name, dim3id, dim3len, &
       dim4name, dim4id, dim4len, &
       dim5name, dim5id, dim5len, &
       dim6name, dim6id, dim6len, &
       dim7name, dim7id, dim7len, &
       dim8name, dim8id, dim8len, &
       dim9name, dim9id, dim9len, &
       dim10name, dim10id, dim10len, &
       dim11name, dim11id, dim11len, &
       dim12name, dim12id, dim12len, &
       dim13name, dim13id, dim13len, &
       Name1,String1, &
       Name2,String2, &
       Name3,String3, &
       Name4,String4, &
       Name5,String5, &
       dbg)
    !
    ! Function that creates a netecf file intended to later fill with
    ! data. Includes global data and saves all the dimensions one
    ! might need. Opens and closes the created netcdf file in
    ! this subroutine

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in) :: file_name        
    INTEGER :: ncid

    !
    ! ------------------------------
    ! FILE_NAME = full path to file
    ! dimXname = the name of a dimension you want to add
    ! dimXid = the dimid that will be associated with the dimension
    ! dimXlen = The len of the dimension


    CHARACTER(len=*), INTENT(in), OPTIONAL :: Name1,Name2,Name3,Name4,Name5
    CHARACTER(len=*), INTENT(in), OPTIONAL :: String1,String2,String3,String4,String5
    CHARACTER(len=*), INTENT(in), OPTIONAL :: dimtname
    CHARACTER(len=*), INTENT(in), OPTIONAL :: dim1name,dim2name,dim3name,dim4name,dim5name
    CHARACTER(len=*), INTENT(in), OPTIONAL :: dim6name,dim7name,dim8name,dim9name,dim10name
    CHARACTER(len=*), INTENT(in), OPTIONAL :: dim11name
    CHARACTER(len=*), INTENT(in), OPTIONAL :: dim12name,dim13name
    INTEGER, INTENT(OUT), OPTIONAL :: dimtid
    INTEGER, INTENT(OUT), OPTIONAL :: dim1id,dim2id,dim3id,dim4id,dim5id
    INTEGER, INTENT(OUT), OPTIONAL :: dim6id,dim7id,dim8id,dim9id,dim10id
    INTEGER, INTENT(OUT), OPTIONAL :: dim11id
    INTEGER, INTENT(OUT), OPTIONAL :: dim12id,dim13id
    INTEGER, INTENT(in), OPTIONAL :: dim1len,dim2len,dim3len,dim4len,dim5len
    INTEGER, INTENT(in), OPTIONAL :: dim6len,dim7len,dim8len,dim9len, dim10len
    INTEGER, INTENT(in), OPTIONAL :: dim11len
    INTEGER, INTENT(in), OPTIONAL :: dim12len,dim13len
    INTEGER, INTENT(in), OPTIONAL :: dbg

    ! Used for saving dimension data
    INTEGER, DIMENSION(8) :: date_time
    INTEGER, DIMENSION(2) :: utc
    CHARACTER(len=1000)   :: create_date

    ncid = 0

    ! Create the file. 
    CALL CHECK(NF90_CREATE(file_name,NF90_CLOBBER,ncid))

    CALL DATE_AND_TIME(values=date_time)

    utc(1) = FLOOR((REAL(60*date_time(5)+date_time(6)-date_time(4)))/60)
    utc(2) = 60*date_time(5)+date_time(6)-date_time(4)  - 60*utc(1)

    WRITE(create_date,'(I4,A1,I2,A1,I2,3x,I2,A1,I2,2x,A3)')&
         date_time(1),'-',date_time(2),'-',date_time(3),&
         utc(1),':',utc(2),'UTC'

    ! write some global attributes
    CALL CHECK(nf90_put_att(ncid,nf90_global,"Contact" ,"Salomon.Eliasson@smhi.se"))
    CALL CHECK(nf90_put_att(ncid,nf90_global,"Creation date",create_date))
    CALL CHECK(nf90_put_att(ncid,nf90_global,"football team","Man U"))
    IF (PRESENT(Name1)) CALL CHECK(NF90_PUT_ATT(ncid,NF90_GLOBAL,TRIM(Name1),TRIM(String1)))
    IF (PRESENT(Name2)) CALL CHECK(NF90_PUT_ATT(ncid,NF90_GLOBAL,TRIM(Name2),TRIM(String2)))
    IF (PRESENT(Name3)) CALL CHECK(NF90_PUT_ATT(ncid,NF90_GLOBAL,TRIM(Name3),TRIM(String3)))
    IF (PRESENT(Name4)) CALL CHECK(NF90_PUT_ATT(ncid,NF90_GLOBAL,TRIM(Name4),TRIM(String4)))
    IF (PRESENT(Name5)) CALL CHECK(NF90_PUT_ATT(ncid,NF90_GLOBAL,TRIM(Name5),TRIM(String5)))

    ! Define the dimensions.
    !
    ! The record dimension is defined to have
    ! unlimited len - it can grow as needed. In this example it is
    ! the time dimension.

    !   define time with unlimited length
    IF (PRESENT(dimtname)) CALL CHECK(NF90_DEF_DIM(ncid,dimtname,NF90_UNLIMITED,dimtid))
    IF (PRESENT(dim1name)) CALL CHECK(NF90_DEF_DIM(ncid,dim1name,dim1len,dim1id))
    IF (PRESENT(dim2name)) CALL CHECK(NF90_DEF_DIM(ncid,dim2name,dim2len,dim2id))
    IF (PRESENT(dim3name)) CALL CHECK(NF90_DEF_DIM(ncid,dim3name,dim3len,dim3id))
    IF (PRESENT(dim4name)) CALL CHECK(NF90_DEF_DIM(ncid,dim4name,dim4len,dim4id))
    IF (PRESENT(dim5name)) CALL CHECK(NF90_DEF_DIM(ncid,dim5name,dim5len,dim5id))
    IF (PRESENT(dim6name).AND.dim6len>0) &
         CALL CHECK(NF90_DEF_DIM(ncid,dim6name,dim6len,dim6id))
    IF (PRESENT(dim7name)) CALL CHECK(NF90_DEF_DIM(ncid,dim7name,dim7len,dim7id))
    IF (PRESENT(dim8name)) CALL CHECK(NF90_DEF_DIM(ncid,dim8name,dim8len,dim8id))
    IF (PRESENT(dim9name)) CALL CHECK(NF90_DEF_DIM(ncid,dim9name,dim9len,dim9id))
    IF (PRESENT(dim10name).AND.dim10len>0) CALL CHECK(NF90_DEF_DIM(ncid,dim10name,dim10len,dim10id))
    IF (PRESENT(dim11name)) CALL CHECK(NF90_DEF_DIM(ncid,dim11name,dim11len,dim11id))
    IF (PRESENT(dim12name)) CALL CHECK(NF90_DEF_DIM(ncid,dim12name,dim12len,dim12id))
    IF (PRESENT(dim13name)) CALL CHECK(NF90_DEF_DIM(ncid,dim13name,dim13len,dim13id))

    CALL CHECK(NF90_ENDDEF(ncid))

    ! Close the file. This causes netCDF to flush all buffers and make
    ! sure your data are really written to disk.
    CALL CHECK(NF90_CLOSE(ncid))


    IF (PRESENT(dbg)) THEN
       IF (dbg > 0 ) PRINT "(A,1x,A)","Created netcdf file with global attributes:",TRIM(file_name)
    END IF
  END SUBROUTINE PREAMBLE_NETCDF

  ! ===============================
  ! Attributes for variables

  SUBROUTINE ADD_VARIABLE_ATTRIBUTES(ncid,varid,name,dimids,ndims,&
       add_offset,coordinates,description,fillvalue_i,fillvalue_r,&
       grid_mapping,isinteger,long_name,scale_factor,unit,valid_min,valid_max,&
       dbg)

    ! Internal function that adds attributes related to variable to be saved to an
    ! already existing and OPEN file created by "preamble_Netcdf".
    !  Also, this function will not close the netcdf file. The main
    !  purpose of this subroutine is to minimise the code repetition
    !  that would be necessary without it

    IMPLICIT NONE

    ! PROVIDE:
    ! file = "%s" fullpath to netcdf file to append to
    ! data = 1 dimensional data
    ! name = "%s" name of short name data field
    ! unit = "%s" the data unit
    ! dimid = "%i" e.g. from nf90_inq_dimid(ncid,"records",dimid)
    ! ndims = number of data dimensions
    !
    ! OPTIONAL:
    ! long_name = "%s" A more human readable name
    ! description = "%s" Description of variable
    ! isinteger = true if your data is integer (but pass it in as real(int_data))

    CHARACTER (len=*), INTENT(in) :: name
    INTEGER, INTENT(in)           :: ncid 
    INTEGER, INTENT(out)          :: varid
    INTEGER, INTENT(in)           :: ndims
    INTEGER, INTENT(in)           :: dimids(ndims)

    CHARACTER (len=*), INTENT(in), OPTIONAL :: long_name,unit,description,coordinates,grid_mapping
    LOGICAL, INTENT(in), OPTIONAL           :: isinteger
    REAL(4), INTENT(in), OPTIONAL           :: fillvalue_r,scale_factor,add_offset,&
                                               valid_min,valid_max
    INTEGER, INTENT(in), OPTIONAL           :: fillvalue_i
    INTEGER, INTENT(in), OPTIONAL           :: dbg

    ! -----> Put the data variable attributes in the file
    IF (PRESENT(isinteger) .AND. isinteger ) THEN
       CALL CHECK(NF90_DEF_VAR(ncid,name,nf90_int,dimids,varid))
    ELSE
       CALL CHECK(NF90_DEF_VAR(ncid,name,nf90_real,dimids,varid))
    END IF
    IF (PRESENT(long_name))     CALL CHECK(NF90_PUT_ATT(ncid,varid,"long_name",long_name))
    IF (PRESENT(description))   CALL CHECK(NF90_PUT_ATT(ncid,varid,"description",description))
    IF (PRESENT(unit))          CALL CHECK(NF90_PUT_ATT(ncid,varid,"units",unit))
    IF (PRESENT(coordinates).AND.TRIM(coordinates).NE.'') &
                                CALL CHECK(NF90_PUT_ATT(ncid,varid,"coordinates",coordinates))
    IF (isinteger) THEN
       IF (PRESENT(fillvalue_i))CALL CHECK(NF90_PUT_ATT(ncid,varid,"_FillValue",fillvalue_i))
    ELSE
       IF (PRESENT(fillvalue_r))CALL CHECK(NF90_PUT_ATT(ncid,varid,"_FillValue",fillvalue_r))
    ENDIF
    IF (PRESENT(grid_mapping).and.TRIM(grid_mapping).NE.'') &
                                CALL CHECK(NF90_PUT_ATT(ncid,varid,"grid_mapping",grid_mapping))
    IF (PRESENT(scale_factor))  CALL CHECK(NF90_PUT_ATT(ncid,varid,"scale_factor",scale_factor))
    IF (PRESENT(add_offset))    CALL CHECK(NF90_PUT_ATT(ncid,varid,"add_offset",add_offset))
    IF (PRESENT(valid_min)) THEN
       IF (valid_min .NE. fillvalue_r) THEN
          CALL CHECK(NF90_PUT_ATT(ncid,varid,"valid_min",valid_min))
       END IF
    END IF
    IF (PRESENT(valid_max)) THEN
       IF (valid_min .NE. fillvalue_r) THEN
          CALL CHECK(NF90_PUT_ATT(ncid,varid,"valid_max",valid_max))
       END IF
    END IF
    CALL CHECK(NF90_ENDDEF(ncid))

    IF (PRESENT(dbg)) THEN
       IF (dbg > 0 ) PRINT '(A,1x,A)',"Added variable attributes for:",name
    END IF

  END SUBROUTINE ADD_VARIABLE_ATTRIBUTES

  ! ---------------------------
  ! check

  SUBROUTINE CHECK(status,dbg,string,fail)

    INTEGER, INTENT (in) :: status
    CHARACTER (len=*), OPTIONAL, INTENT(in) :: string
    INTEGER, OPTIONAL, INTENT(in):: dbg
    LOGICAL, OPTIONAL, INTENT(out) :: fail
    LOGICAL :: hasdbg,hasfail
    
    hasdbg = PRESENT(dbg) .AND. PRESENT(string)
    hasfail=.FALSE. ! reusing this
    IF ((hasdbg) .AND. (dbg>0)) THEN
       PRINT *,string
    END IF

    IF (status /= nf90_noerr) THEN
       PRINT *,TRIM(nf90_strerror(status))
       hasfail=.TRUE.
    END IF
    IF (PRESENT(fail)) fail=hasfail

  END SUBROUTINE CHECK

END MODULE MY_NETCDFTOOLS
