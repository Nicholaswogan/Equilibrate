!> ------------------------------------
!>        EASY CHEM, a CEA clone
!> ------------------------------------


!> DATA STORAGE MODULE
module data_block
   implicit none
   ! private

   ! public :: SET_DATA
   ! public :: EASYCHEM

   !> Run variables
   integer     :: iter_max
   integer     :: N_atoms, N_reactants, N_gas, N_cond, N_ions
   logical     :: verbose, verbose_cond, quick, ions, remove_ions, error
   character(len=500)   :: err_msg

   !> Constants
   double precision, parameter :: R = 8.3144598d0
   double precision, parameter :: amu = 1.660538921d-24
   double precision, parameter :: kB = 1.3806488d-16
   double precision, parameter :: mol = 6.02214129d23

   !> List of atoms
   character(len=2), allocatable   :: names_atoms(:)  !(N_atoms)
   integer, allocatable            :: id_atoms(:)  !(N_atoms)

   !> List of reactants reordered with condensates at the end
   character(len=15), allocatable  :: names_reactants(:), names_reactants_orig(:)  !(N_reac)
   integer, allocatable            :: id_reactants(:,:)  !(N_reac,2)

   !> Atomic data for each reactant
   character(len=2), allocatable   :: reac_atoms_names(:,:)  !(5,N_reac)
   integer, allocatable            :: reac_atoms_id(:,:)  !(5,N_reac)
   double precision, allocatable   :: reac_stoich(:,:)  !(5,N_reac)

   !> Nature of reactant
   logical, allocatable    :: reac_condensed(:), reac_ion(:)  !(N_reac)

   !> Thermodynamic data arrays
   integer, parameter      :: N_coeffs = 10, N_temps = 10
   integer, allocatable    :: thermo_data_n_coeffs(:,:)  !(N_temps,N_reac)
   integer, allocatable    :: thermo_data_n_intervs(:)  !(N_reac)
   double precision, allocatable   :: thermo_data(:,:,:)  !(N_coeffs,N_temps,N_reac)
   double precision, allocatable   :: thermo_data_temps(:,:,:)  !(2,N_temps,N_reac)
   double precision, allocatable   :: thermo_data_T_exps(:,:,:)  !(8,N_temps,N_reac)
   double precision, allocatable   :: form_heat_Jmol_298_15_K(:)  !(N_reac)
   double precision, allocatable   :: H_0_298_15_K_m_H_0_0_K(:,:)  !(N_temps, N_reac)
   double precision, allocatable   :: mol_weight(:)  !(N_reac)

   !> Atoms & masses, from http://www.science.co.il/PTelements.asp
   integer, parameter      :: N_atoms_save = 104
   character*2, parameter  :: names_atoms_save(N_atoms_save) = &
   (/ 'E ','H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na', &
   'Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn', &
   'Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr', &
   'Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb', &
   'Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd', &
   'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir', &
   'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
   'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr' /)
   double precision, parameter  :: masses_atoms_save(N_atoms_save) = &
   amu*(/ 0.000548579909d0, 1.0079d0,4.0026d0,6.941d0,9.0122d0,10.811d0,12.0107d0,14.0067d0 &
   ,15.9994d0,18.9984d0,20.1797d0,22.9897d0,24.305d0,26.9815d0,28.0855d0,30.9738d0,32.065d0 &
   ,35.453d0,39.948d0,39.0983d0,40.078d0,44.9559d0,47.867d0,50.9415d0,51.9961d0,54.938d0,55.845d0 &
   ,58.9332d0,58.6934d0,63.546d0,65.39d0,69.723d0,72.64d0,74.9216d0,78.96d0,79.904d0,83.8d0,85.4678d0 &
   ,87.62d0,88.9059d0,91.224d0,92.9064d0,95.94d0,98d0,101.07d0,102.9055d0,106.42d0,107.8682d0 &
   ,112.411d0,114.818d0,118.71d0,121.76d0,127.6d0,126.9045d0,131.293d0,132.9055d0,137.327d0,138.9055d0 &
   ,140.116d0,140.9077d0,144.24d0,145d0,150.36d0,151.964d0,157.25d0,158.9253d0,162.5d0,164.9303d0 &
   ,167.259d0,168.9342d0,173.04d0,174.967d0,178.49d0,180.9479d0,183.84d0,186.207d0,190.23d0,192.217d0 &
   ,195.078d0,196.9665d0,200.59d0,204.3833d0,207.2d0,208.9804d0,209d0,210d0,222d0,223d0,226d0,227d0,232.0381d0 &
   ,231.0359d0,238.0289d0,237d0,244d0,243d0,247d0,247d0,251d0,252d0,257d0,258d0,259d0,262d0/)

   contains

   !> INITIALIZE ALL DATA
   subroutine SET_DATA(N_atoms_in, N_reactants_in, atoms_char, reac_char, fpath)
      character, intent(in)   :: atoms_char(N_atoms_in,2), reac_char(N_reactants_in,15)
      integer, intent(in)     :: N_atoms_in, N_reactants_in
      character(len=800), intent(in) :: fpath

      interface
      subroutine chrarr_to_stringarr(chr, str)
         character, intent(in)                   :: chr(:,:)         ! row = 1 string
         character(len=size(chr,2)), intent(out) :: str(size(chr,1))
      end subroutine
      end interface

      error = .false.

      ! REACTANTS
      if (N_reactants_in /= N_reactants .and. allocated(names_reactants_orig)) then
         ! Deallocate everything if the number of reactants changed
         deallocate(names_reactants_orig)
         deallocate(names_reactants)
         deallocate(id_reactants)
         deallocate(reac_atoms_names)
         deallocate(reac_atoms_id)
         deallocate(reac_stoich)
         deallocate(reac_condensed)
         deallocate(reac_ion)
         deallocate(thermo_data_n_coeffs)
         deallocate(thermo_data_n_intervs)
         deallocate(thermo_data)
         deallocate(thermo_data_temps)
         deallocate(thermo_data_T_exps)
         deallocate(form_heat_Jmol_298_15_K)
         deallocate(H_0_298_15_K_m_H_0_0_K)
         deallocate(mol_weight)
      end if
      if (.not. allocated(names_reactants_orig)) then
         ! Allocate reactant-related arrays
         N_reactants = N_reactants_in
         allocate(names_reactants_orig(N_reactants))
         allocate(names_reactants(N_reactants))
         allocate(id_reactants(N_reactants,2))
         allocate(reac_atoms_names(5,N_reactants))
         allocate(reac_atoms_id(5,N_reactants))
         allocate(reac_stoich(5,N_reactants))
         allocate(reac_condensed(N_reactants))
         allocate(reac_ion(N_reactants))

         allocate(thermo_data_n_coeffs(N_temps,N_reactants))
         allocate(thermo_data_n_intervs(N_reactants))
         allocate(thermo_data(N_coeffs,N_temps,N_reactants))
         allocate(thermo_data_temps(2,N_temps,N_reactants))
         allocate(thermo_data_T_exps(8,N_temps,N_reactants))
         allocate(form_heat_Jmol_298_15_K(N_reactants))
         allocate(H_0_298_15_K_m_H_0_0_K(N_temps,N_reactants))
         allocate(mol_weight(N_reactants))
      end if

      ! Set names_reactants_orig with the given list of reactants
      call da_CH2STR(reac_char, names_reactants_orig)

      ! Set all thermodynamic data
      call da_READ_THERMO(fpath)
      if (error) RETURN

      ! ATOMS
      if (N_atoms_in+1 /= N_atoms .and. allocated(names_atoms)) then
         deallocate(names_atoms)
         deallocate(id_atoms)
      end if
      if (.not. allocated(names_atoms)) then
         N_atoms = N_atoms_in + 1
         allocate(names_atoms(N_atoms))
         allocate(id_atoms(N_atoms))
      end if

      call da_CH2STR(atoms_char, names_atoms(1:N_atoms_in))
      names_atoms(N_atoms) = 'E'

      call da_ATOMS_ID()
      if (error) RETURN

      call da_REAC_ATOMS_ID()

   end subroutine SET_DATA

   !> Convert a 2-dimensional array of characters into an array of strings ; used in SET_DATA
   subroutine da_CH2STR(chr, str)

      character, intent(in)                   :: chr(:,:)  ! 1 row = 1 string
      character(len=size(chr,2)), intent(out) :: str(size(chr,1))
      integer                                 :: i, j

      str = ''
      do i = 1, size(chr,1)
         do j = 1, size(chr,2)
            str(i) = trim(str(i))//chr(i,j)
         end do
      end do
   end subroutine da_CH2STR

   !> Sets id_atoms, where the i-th cell corresponds to names_atoms(i) and contains the index of the same atom in names_atoms_save
   !> names_atoms_save(id_atoms(i)) = names_atoms(i)
   subroutine da_ATOMS_ID()

      integer           :: i_atom, i_atom_save
      logical           :: change
      character(len=2)  :: atom_upper, atom_upper_save
      character(len=3)  :: num

      do i_atom = 1, N_atoms
         if (trim(names_atoms(i_atom)) == '') then
            id_atoms(i_atom) = -1
            CYCLE
         end if
         change = .false.
         call uppercase(names_atoms(i_atom), atom_upper)
         do i_atom_save = 1, N_atoms_save
            call uppercase(names_atoms_save(i_atom_save), atom_upper_save)
            if (atom_upper == atom_upper_save) then
               id_atoms(i_atom) = i_atom_save
               change = .true.
               exit
            end if
         end do
         if (.not. change) then
            id_atoms(i_atom) = -1
            error = .true.
            write(num, '(i3.3)') i_atom
            err_msg = "READ DATA ERROR: the atom #"//num//" '"//names_atoms(i_atom)//"' was not recognised."
            RETURN
         end if
      end do
   end subroutine da_ATOMS_ID

   subroutine da_REAC_ATOMS_ID()

      integer           :: i_reac, i_atom, i_atom_save
      logical           :: change
      character(len=2)  :: atom_upper, atom_upper_save
      character(len=3)  :: num

      do i_reac = 1, N_reactants
         do i_atom = 1, 5
            if (trim(reac_atoms_names(i_atom,i_reac))=='') then
               reac_atoms_id(i_atom, i_reac) = -1
               CYCLE
            end if
            change = .false.
            call uppercase(reac_atoms_names(i_atom,i_reac), atom_upper)
            do i_atom_save = 1, N_atoms_save
               call uppercase(names_atoms_save(i_atom_save), atom_upper_save)
               if (atom_upper == atom_upper_save) then
                  reac_atoms_id(i_atom,i_reac) = i_atom_save
                  change = .true.
                  EXIT
               end if
            end do
            if (.not. change) then
               reac_atoms_id(i_atom, i_reac) = -1
               error = .true.
               write(num, '(i3.3)') i_atom
               err_msg = "READ DATA ERROR: the atom #"//num//" '"//names_atoms(i_atom)//&
               "' in reactant '"//names_reactants(i_reac)//"' was not recognised."
               RETURN
            end if
         end do
      end do
   end subroutine da_REAC_ATOMS_ID

   !> Read in provided file all thermodynamic data
   subroutine da_READ_THERMO(fpath)

      character(len=800), intent(in) :: fpath
      character(len=80)             :: file_line, file_line_up
      character(len=15)             :: name_reac_up
      integer                       :: i_reac, n_interv, i_interv, i_stoich, stoich_start, reac_found

      reac_found = 0
      thermo_data_n_intervs = -1
      N_gas = 0
      N_ions = 0

      reac_ion = .FALSE.
      ions = .FALSE.

      open(unit=17,file=fpath, action="READ", status="OLD")
      do while (1>0)
         read(17,'(A80)',end=122) file_line
         call uppercase(file_line, file_line_up)
         do i_reac = 1, N_reactants
            call uppercase(names_reactants_orig(i_reac), name_reac_up)

            if (trim(adjustl(file_line_up(1:18))) == trim(adjustl(name_reac_up))) then
               ! Recognized a reactant in file
               reac_found = reac_found+1

               ! Mark it as found
               read(17,'(A80)',end=122) file_line
               read(file_line(1:3),'(I2)') n_interv
               thermo_data_n_intervs(i_reac) = n_interv

               ! Gas or condensate ?
               if (file_line(52:52) == '0') then
                  reac_condensed(i_reac) = .FALSE.
                  N_gas = N_gas + 1
               else
                  reac_condensed(i_reac) = .TRUE.
               end if
            end if

         end do
      end do
      122 close(17)

      ! If not all reactants were found
      if (reac_found /= N_reactants) then
         error = .true.
         err_msg = 'READ DATA ERROR: For the following species no thermodynamical data was found:'
         n_interv = 0
         do i_reac = 1, N_reactants
            if (thermo_data_n_intervs(i_reac) .EQ. -1) then
               if (n_interv == 0) then
                  err_msg = trim(err_msg)//' '//trim(adjustl(names_reactants_orig(i_reac)))
               else
                  err_msg = trim(err_msg)//', '//trim(adjustl(names_reactants_orig(i_reac)))
               end if
            end if
         end do
         RETURN
      end if

      N_cond = N_reactants - N_gas

      ! Puts in names_reactants the ordered list of reactants, sets id_reactants as a link between the two
      call da_REORDER_SPECS()

      ! BASED ON THE ~RIGHT DESCRIPTION GIVEN IN GORDON 1996, page 73 AND THE APPEARANCE OF THERMO.INP.
      open(unit=17,file=fpath, action="READ", status="OLD")
      do while (1>0)
         read(17,'(A80)',end=123) file_line
         call uppercase(file_line, file_line_up)
         do i_reac = 1, N_reactants
            call uppercase(names_reactants(i_reac), name_reac_up)

            if (trim(adjustl(file_line_up(1:18))) == trim(adjustl(name_reac_up))) then
               ! Recognized a reactant in file
               read(17,'(A80)',end=123) file_line
               read(file_line(1:3),'(I2)') n_interv
               thermo_data_n_intervs(i_reac) = n_interv

               stoich_start = 11
               do i_stoich = 1, 5
                  ! Set the atoms for each reactant
                  reac_atoms_names(i_stoich,i_reac) = file_line(stoich_start:stoich_start+1)
                  ! Are there ions to be treated?
                  if (trim(adjustl(reac_atoms_names(i_stoich,i_reac))) == 'E') then
                     ions = .TRUE.
                     reac_ion(i_reac) = .TRUE.
                     N_ions = N_ions + 1
                  end if
                  read(file_line(stoich_start+2:stoich_start+7),'(F6.2)') reac_stoich(i_stoich,i_reac)
                  stoich_start = stoich_start+8
               end do

               read(file_line(53:65),'(F13.5)') mol_weight(i_reac)
               read(file_line(66:80),'(F13.5)') form_heat_Jmol_298_15_K(i_reac)

               do i_interv = 1, 3*n_interv
                  read(17,'(A80)',end=123) file_line
                  if (MOD(i_interv,3) == 1) then
                     read(file_line(2:22),'(F10.3,X,F10.3)') thermo_data_temps(1,i_interv/3+1,i_reac), &
                     thermo_data_temps(2,i_interv/3+1,i_reac)
                     read(file_line(23:23),'(I1)') thermo_data_n_coeffs(i_interv/3+1,i_reac)
                     read(file_line(24:63),'(8F5.1)') thermo_data_T_exps(1,i_interv/3+1,i_reac), &
                     thermo_data_T_exps(2,i_interv/3+1,i_reac), thermo_data_T_exps(3,i_interv/3+1,i_reac), &
                     thermo_data_T_exps(4,i_interv/3+1,i_reac), thermo_data_T_exps(5,i_interv/3+1,i_reac), &
                     thermo_data_T_exps(6,i_interv/3+1,i_reac), thermo_data_T_exps(7,i_interv/3+1,i_reac), &
                     thermo_data_T_exps(8,i_interv/3+1,i_reac)
                     read(file_line(66:80),'(F15.3)') H_0_298_15_K_m_H_0_0_K(i_interv/3+1,i_reac)
                  end if
                  if (MOD(i_interv,3) == 2) then
                     read(file_line(1:80),'(5D16.8)') thermo_data(1,i_interv/3+1,i_reac), &
                     thermo_data(2,i_interv/3+1,i_reac), thermo_data(3,i_interv/3+1,i_reac) , &
                     thermo_data(4,i_interv/3+1,i_reac), thermo_data(5,i_interv/3+1,i_reac)
                  end if
                  if (MOD(i_interv,3) == 0) then
                     read(file_line(1:80),'(5D16.8)') thermo_data(6,i_interv/3,i_reac), &
                     thermo_data(7,i_interv/3,i_reac), thermo_data(8,i_interv/3,i_reac) , &
                     thermo_data(9,i_interv/3,i_reac), thermo_data(10,i_interv/3,i_reac)
                  end if
               end do
            end if
         end do
      end do
      123 CLOSE(17)
   end subroutine da_READ_THERMO

   !> Puts every letter in strin in uppercase
   subroutine uppercase(strin,strout)

      character(len=*)  :: strin,strout
      integer           :: i

      strout = strin
      do i = 1, len(strout)
         select case(strout(i:i))
         case("a":"z")
            strout(i:i) = achar(iachar(strout(i:i))-32)
         end select
      end do
   end subroutine uppercase

   !> Puts in names_reactants the ordered list of reactants, sets id_reactants as a link between the two
   subroutine da_REORDER_SPECS()

      integer  :: i_reac, gas_offset, cond_offset

      gas_offset = 1
      cond_offset = 1
      do i_reac = 1, N_reactants
         if (reac_condensed(i_reac)) then
            names_reactants(N_gas+cond_offset) = names_reactants_orig(i_reac)
            ! 1st row = new index of the original reactant at i_reac
            id_reactants(i_reac,1) = N_gas + cond_offset
            ! 2nd row = original index of reactants in ordered list
            id_reactants(N_gas+cond_offset,2) = i_reac
            cond_offset = cond_offset + 1
         else
            names_reactants(gas_offset) = names_reactants_orig(i_reac)
            id_reactants(i_reac,1) = gas_offset
            id_reactants(gas_offset,2) = i_reac
            gas_offset = gas_offset + 1
         end if
      end do
   end subroutine da_REORDER_SPECS



   !> MAIN SUBROUTINE
   subroutine EASYCHEM(mode,verbo,N_atoms_in,N_reactants_in,molfracs_atoms, &
      molfracs_reactants,massfracs_reactants,temp,press,nabla_ad,gamma2,MMW,rho,c_pe)

      !! I/O:
      character, intent(in)            :: mode
      character(len=2), intent(in)     :: verbo
      double precision, intent(in)     :: molfracs_atoms(N_atoms_in)
      double precision, intent(out)    :: molfracs_reactants(N_reactants_in), massfracs_reactants(N_reactants_in)
      integer, intent(in)              :: N_atoms_in, N_reactants_in
      double precision, intent(in)     :: temp, press
      double precision, intent(out)    :: nabla_ad,gamma2,MMW,rho,c_pe

      !! Internal:
      double precision                 :: C_P_0(N_reactants), H_0(N_reactants), S_0(N_reactants)
      double precision                 :: molfracs_atoms_ions(N_atoms_in+1), temp_use
      integer                          :: N_atoms_use, gamma_neg_try

      error = .false.

      if (N_atoms /= N_atoms_in+1 .or. N_reactants /= N_reactants_in) then
         error = .true.
         err_msg = "VALUE ERROR: The initialized and given arrays are not of the same size..."
         RETURN
      end if

      verbose = .FALSE.
      verbose_cond = (verbo == 'vy')
      quick = (mode /= 's')
      remove_ions = .FALSE.

      call INIT_RAND_SEED()

      molfracs_atoms_ions(1:N_atoms_in) = molfracs_atoms
      molfracs_atoms_ions(N_atoms) = 0d0

      if (ions) then
         if (temp > 750) then
            N_atoms_use = N_atoms
         else
            remove_ions = .true.
            N_atoms_use = N_atoms_in
         end if
      else
         N_atoms_use = N_atoms_in
      end if

      ! CALCULATION BEGINS

      call ec_COMP_THERMO_QUANTS(temp,N_reactants,C_P_0, H_0, S_0)
      gamma2 = 0d0
      temp_use = temp
      gamma_neg_try = 0d0
      do while (gamma2 < 1d0)
         call ec_COMP_EQU_CHEM(N_atoms_use, N_reactants,  molfracs_atoms_ions(1:N_atoms_use), &
         molfracs_reactants, massfracs_reactants, &
         temp_use, press, C_P_0, H_0, S_0, &
         nabla_ad, gamma2, MMW, rho, c_pe)
         if (error) RETURN

         if (gamma2 < 1d0) then
            write(*,*) 'Gamma was < 1, redo! gamma2, temp, ', gamma2, temp
            gamma_neg_try = gamma_neg_try + 1
            if (gamma_neg_try > 10) then
               call random_number(temp_use)
               temp_use = temp*(1d0 + 0.01d0*temp_use)
               write(*,*) 'temp, temp_use', temp, temp_use
               call ec_COMP_THERMO_QUANTS(temp_use,N_reactants,C_P_0, H_0, S_0)
            end if
         end if
      end do

      c_pe = c_pe*1d7 ! J/(g K) to erg/(g K)

   end subroutine EASYCHEM

   !> Computes the values of C_P_0, H_0 and S_0
   subroutine ec_COMP_THERMO_QUANTS(temp,N_reac,C_P_0, H_0, S_0)
      !! I/O
      double precision, intent(in)  :: temp
      integer, intent(in)           :: N_reac
      double precision, intent(out) :: C_P_0(N_reac), H_0(N_reac), &
      S_0(N_reac)
      !! internal
      integer                       :: i_reac, i_tempinv, tempinv_ind

      do i_reac = 1, N_reactants

         ! Get temperature interpolation range
         if (temp < thermo_data_temps(1,1,i_reac)) then
            tempinv_ind = 1
         else if (temp >= thermo_data_temps(2,thermo_data_n_intervs(i_reac),i_reac)) then
            tempinv_ind = thermo_data_n_intervs(i_reac)
         else
            do i_tempinv = 1, thermo_data_n_intervs(i_reac)
               if ((temp >= thermo_data_temps(1,i_tempinv,i_reac)) .AND. &
               (temp < thermo_data_temps(2,i_tempinv,i_reac))) then
                  tempinv_ind = i_tempinv
                  EXIT
               end if
            end do
         end if

         ! Calculate thermodynamic quantities as explained in Gordon 1996, page 74
         C_P_0(i_reac) = (thermo_data(1,tempinv_ind,i_reac)*temp**(-2)+ &
         thermo_data(2,tempinv_ind,i_reac)*temp**(-1)+ &
         thermo_data(3,tempinv_ind,i_reac)+thermo_data(4,tempinv_ind,i_reac)* &
         temp**(1)+thermo_data(5,tempinv_ind,i_reac)*temp**(2)+ &
         thermo_data(6,tempinv_ind,i_reac)*temp**(3)+thermo_data(7,tempinv_ind,i_reac)* &
         temp**(4))*R
         H_0(i_reac) = (-thermo_data(1,tempinv_ind,i_reac)*temp**(-2)+ &
         thermo_data(2,tempinv_ind,i_reac)*temp**(-1)*log(temp)+ &
         thermo_data(3,tempinv_ind,i_reac)+thermo_data(4,tempinv_ind,i_reac)*temp**(1)/2d0+ &
         thermo_data(5,tempinv_ind,i_reac)*temp**(2)/3d0+ &
         thermo_data(6,tempinv_ind,i_reac)*temp**(3)/4d0+thermo_data(7,tempinv_ind,i_reac)* &
         temp**(4)/5d0+thermo_data(9,tempinv_ind,i_reac)/temp)* &
         R*temp
         S_0(i_reac) = (-thermo_data(1,tempinv_ind,i_reac)*temp**(-2)/2d0- &
         thermo_data(2,tempinv_ind,i_reac)*temp**(-1)+ &
         thermo_data(3,tempinv_ind,i_reac)*log(temp)+ &
         thermo_data(4,tempinv_ind,i_reac)*temp**(1)+ &
         thermo_data(5,tempinv_ind,i_reac)*temp**(2)/2d0+ &
         thermo_data(6,tempinv_ind,i_reac)*temp**(3)/3d0+thermo_data(7,tempinv_ind,i_reac)* &
         temp**(4)/4d0+thermo_data(10,tempinv_ind,i_reac))*R

      end do

   end subroutine ec_COMP_THERMO_QUANTS

   !> Computes the specie abundances (molar and mass)
   recursive subroutine ec_COMP_EQU_CHEM(N_atoms_use, N_reac, molfracs_atoms, &
      molfracs_reactants, massfracs_reactants, &
      temp, press, C_P_0, H_0, S_0, nabla_ad, gamma2, MMW, rho, c_pe)

      !! I/O:
      integer, intent(in)              :: N_atoms_use, N_reac
      double precision, intent(in)     :: molfracs_atoms(N_atoms_use)
      double precision, intent(inout)  :: molfracs_reactants(N_reac), massfracs_reactants(N_reac)
      DOUBLE PRECISION, intent(in)     :: temp, press
      DOUBLE PRECISION, intent(in)     :: C_P_0(N_reac), H_0(N_reac), S_0(N_reac)
      double precision, intent(out)    :: nabla_ad, gamma2, MMW, rho, c_pe

      !! CEA McBride 1994 style variables:
      DOUBLE PRECISION  :: n ! Moles of gas particles per total mass of mixture in kg
      DOUBLE PRECISION  :: n_spec(N_reactants) ! Moles of species per total mass of mixture in kg
      DOUBLE PRECISION  :: n_spec_old(N_reactants) ! Moles of species per total mass of mixture in kg of previous iteration
      DOUBLE PRECISION  :: pi_atom(N_atoms_use) ! Lagrangian multipliers for the atomic species divided by (R*T)
      DOUBLE PRECISION  :: matrix(N_reactants+N_atoms_use+1,N_reactants+N_atoms_use+1)
      ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for condensed species, the pis and the delta log(n)
      DOUBLE PRECISION  :: vector(N_reactants+N_atoms_use+1), solution_vector(N_reactants+N_atoms_use+1)

      !! Internal:
      INTEGER           :: i_iter, i_reac, inc_next, current_solids_number, N_spec_eff, buffer_ind, i_atom
      LOGICAL           :: converged, remove_cond, slowed
      LOGICAL           :: solid_inclu(N_cond), neg_cond(N_cond)
      DOUBLE PRECISION  :: dgdnj(N_cond)
      INTEGER           :: solid_indices(N_cond), solid_indices_buff(N_cond)
      DOUBLE PRECISION  :: nsum, mu_gas(N_gas), a_gas(N_gas,N_atoms_use), mass_species, atom_mass, msum

      converged = .FALSE.
      slowed = .FALSE.
      call ec_INIT_ALL_VALS(N_atoms_use,N_reactants,n,n_spec,pi_atom)

      iter_max = 50000 + N_reactants/2
      current_solids_number = 0

      MMW = 0d0

      n_spec_old = n_spec

      ! FIRST: DO GAS ONLY!
      DO i_iter = 1, iter_max

         ! IF (quick) THEN
         call ec_PREP_MATRIX_SHORT(N_atoms_use,N_reactants, molfracs_atoms,N_gas,press,temp,&
         H_0,S_0,n,n_spec,matrix(1:N_atoms_use+1,1:N_atoms_use+1),vector(1:N_atoms_use+1),&
         (/1,1,1,1,1/),5, mu_gas,a_gas)
         call ec_INVERT_MATRIX_SHORT(N_atoms_use+1, &
         matrix(1:N_atoms_use+1,1:N_atoms_use+1),vector(1:N_atoms_use+1), &
         solution_vector(1:N_atoms_use+1))
         if (error) RETURN

         call ec_UPDATE_ABUNDS_SHORT(N_atoms_use,N_reactants,N_gas,&
         solution_vector(1:N_atoms_use+1), n_spec,pi_atom,n,converged,&
         (/1,1,1,1,1/),5,mu_gas,a_gas,temp,molfracs_atoms,n_spec_old)
         ! ELSE
         !    call ec_PREP_MATRIX_LONG(N_atoms,id_atoms,molfracs_atoms,N_gas,&
         !    press,temp,H_0,S_0,n,n_spec,&
         !    matrix(1:N_gas+N_atoms+1,1:N_gas+N_atoms+1),vector(1:N_gas+N_atoms+1),&
         !    (/1,1,1,1,1/),names_reactants,N_reactants,5)
         !    call ec_INVERT_MATRIX_LONG(N_atoms+N_gas+1, &
         !    matrix(1:N_gas+N_atoms+1,1:N_gas+N_atoms+1),vector(1:N_gas+N_atoms+1), &
         !    solution_vector(1:N_gas+N_atoms+1))
         !    call ec_UPDATE_ABUNDS_LONG(N_atoms,N_gas,solution_vector(1:N_gas+N_atoms+1), &
         !    n_spec,pi_atom,n,converged,(/1,1,1,1,1/),5,id_atoms,molfracs_atoms,N_reactants,&
         !    n_spec_old)
         ! END IF

         n_spec_old = n_spec

         IF (verbose) THEN
            write(*,*)
            write(*,*)
            write(*,*) i_iter
            DO i_reac = 1, N_reactants
               write(*,*) names_reactants(i_reac), n_spec(i_reac)/SUM(n_spec)
            END DO
         END IF

         IF (converged) THEN
            EXIT
         END IF
      END DO

      IF (.NOT. converged) THEN
         WRITE(*,*) 'EASY CHEM WARNING: One or more convergence criteria not satisfied! Press, temp', press, temp
         print *
      END IF

      converged = .FALSE.
      remove_cond = .FALSE.

      IF (N_gas .EQ. N_reactants) THEN
         call ec_COMP_ADIABATIC_GRAD(N_atoms_use, N_reactants, N_gas,  n_spec, &
         n, H_0, C_P_0, (/ 1,1,1,1,1 /), 5, temp, nabla_ad, gamma2, c_pe)
         if (error) RETURN
      END IF

      ! THEN: INCLUDE CONDENSATES!
      IF (N_cond > 0) THEN
         solid_inclu = .FALSE.
         inc_next = 0
         neg_cond = .FALSE.

         N_spec_eff = N_gas

         DO WHILE (inc_next /= -1)

            call ec_INCLUDE_WHICH_SOLID(N_atoms_use,N_reac,pi_atom,H_0,S_0,temp, &
            n_spec,solid_inclu,neg_cond,dgdnj,remove_cond,inc_next)

            if (inc_next/=-1) then

               IF (remove_cond) THEN
                  current_solids_number = current_solids_number - 1
                  solid_indices_buff = 0
                  buffer_ind = 1
                  DO i_reac = 1, N_reactants-N_gas
                     IF (solid_indices(i_reac) .NE. inc_next) THEN
                        solid_indices_buff(buffer_ind) = solid_indices(i_reac)
                        buffer_ind = buffer_ind + 1
                     END IF
                  END DO
                  solid_indices = solid_indices_buff
                  solid_inclu(inc_next-N_gas) = .FALSE.
                  neg_cond(inc_next-N_gas) = .TRUE.
                  if (verbose_cond) then
                     print *, '-   ', names_reactants(inc_next), dgdnj(inc_next-N_gas), n_spec(inc_next)
                  end if
                  n_spec(inc_next) = 0d0
               ELSE
                  current_solids_number = current_solids_number + 1
                  solid_indices(current_solids_number) = inc_next
                  solid_inclu(inc_next-N_gas) = .TRUE.
                  call ec_INIT_COND_VALS(N_atoms_use, N_reactants, molfracs_atoms, inc_next, n_spec)
                  if (verbose_cond) then
                     print *, '+   ', names_reactants(inc_next), dgdnj(inc_next-N_gas), n_spec(inc_next)
                  end if
               END IF

               N_spec_eff = N_gas+current_solids_number
               DO i_iter = 1, iter_max

                  IF (quick) THEN
                     call ec_PREP_MATRIX_SHORT(N_atoms_use, N_reactants, molfracs_atoms,N_spec_eff, &
                     press, temp, H_0, S_0, n, n_spec, &
                     matrix(1:N_atoms_use+1+N_spec_eff-N_gas,1:N_atoms_use+1+N_spec_eff-N_gas), &
                     vector(1:N_atoms_use+1+N_spec_eff-N_gas), solid_indices, &
                     N_spec_eff-N_gas, mu_gas, a_gas)
                     call ec_INVERT_MATRIX_SHORT(N_atoms_use+1+N_spec_eff-N_gas, &
                     matrix(1:N_atoms_use+1+N_spec_eff-N_gas,1:N_atoms_use+1+N_spec_eff-N_gas), &
                     vector(1:N_atoms_use+1+N_spec_eff-N_gas), &
                     solution_vector(1:N_atoms_use+1+N_spec_eff-N_gas))
                     if (error) RETURN

                     call ec_UPDATE_ABUNDS_SHORT(N_atoms_use,N_reactants,N_spec_eff,&
                     solution_vector(1:N_atoms_use+1+N_spec_eff-N_gas), &
                     n_spec,pi_atom,n,converged,&
                     solid_indices,N_spec_eff-N_gas,mu_gas,a_gas,temp,molfracs_atoms, &
                     n_spec_old)
                  ELSE
                     call ec_PREP_MATRIX_LONG(N_atoms_use,N_reactants,molfracs_atoms,N_spec_eff,&
                     press,temp,H_0,S_0,n,n_spec,&
                     matrix(1:N_spec_eff+N_atoms_use+1,1:N_spec_eff+N_atoms_use+1),&
                     vector(1:N_spec_eff+N_atoms_use+1),&
                     solid_indices,N_spec_eff-N_gas)
                     call ec_INVERT_MATRIX_LONG(N_atoms_use+N_spec_eff+1, &
                     matrix(1:N_spec_eff+N_atoms_use+1,1:N_spec_eff+N_atoms_use+1),&
                     vector(1:N_spec_eff+N_atoms_use+1), &
                     solution_vector(1:N_spec_eff+N_atoms_use+1))
                     if (error) RETURN

                     call ec_UPDATE_ABUNDS_LONG(N_atoms_use, N_reac, N_spec_eff, &
                     solution_vector(1:N_spec_eff+N_atoms_use+1), &
                     n_spec, pi_atom, n, converged, &
                     solid_indices, N_spec_eff-N_gas, molfracs_atoms, n_spec_old)
                  END IF

                  ! call writetxtall(N_reactants, n_spec)
                  n_spec_old = n_spec

                  IF (verbose) THEN
                     write(*,*)
                     write(*,*)
                     write(*,*) i_iter
                     DO i_reac = 1, N_reactants
                        write(*,*) names_reactants(i_reac), n_spec(i_reac)/SUM(n_spec)
                     END DO
                  END IF

                  DO i_reac = N_gas+1, N_reactants
                     IF ((n_spec(i_reac) < 0d0) .AND. (i_iter > 30)) THEN
                        converged = .TRUE.
                        EXIT
                     END IF
                  END DO

                  IF (converged) THEN
                     EXIT
                  END IF

               END DO

               IF (.NOT. converged) THEN
                  IF (quick) THEN
                     quick = .FALSE.
                     print *
                     print *, 'SLOW ! Press, Temp', press, temp
                     print *
                     call ec_COMP_EQU_CHEM(N_atoms_use, N_reactants, molfracs_atoms, &
                     molfracs_reactants, massfracs_reactants, &
                     temp, press, C_P_0, H_0, S_0, &
                     nabla_ad,gamma2,MMW,rho,c_pe)
                     quick = .TRUE.
                     slowed = .TRUE.
                     EXIT
                  ELSE
                     WRITE(*,*) 'EASY CHEM WARNING: One or more convergence criteria' // &
                                             'not satisfied! in cond Press, temp', press, temp
                     print *
                  END IF
               END IF

               converged = .FALSE.

            END IF

            remove_cond = .FALSE.

         END DO

         ! Calc. nabla_ad
         IF (.NOT. slowed) THEN
            call ec_COMP_ADIABATIC_GRAD(N_atoms_use, N_reactants, N_spec_eff, n_spec, &
            n,H_0,C_P_0,solid_indices,N_spec_eff-N_gas,temp, nabla_ad,gamma2,c_pe)
            if (error) RETURN
         END IF

         if (verbose_cond) then
            print *
            print *, 'Solids included:'
            do i_reac = 1, N_reactants-N_gas
               if (solid_inclu(i_reac)) then
                  print *, ' ', names_reactants(i_reac+N_gas), n_spec(i_reac+N_gas)
               end if
            end do
         end if

      END IF

      ! PREPARE FINAL OUTPUT

      IF (.NOT. slowed) THEN
         nsum = SUM(n_spec)
         DO i_reac = 1, N_reactants
            IF (n_spec(i_reac)/nsum < 1d-50) THEN
               n_spec(i_reac) = 0d0
            END IF
         END DO

         nsum = SUM(n_spec)
         do i_reac = 1, N_reactants
            molfracs_reactants(i_reac) = n_spec(id_reactants(i_reac,1))/nsum
         end do

         msum = 0d0
         DO i_reac = 1, N_reactants
            mass_species = 0d0
            DO i_atom = 1, 5
               if (reac_atoms_id(i_atom,i_reac)>0) then
                  atom_mass = masses_atoms_save(reac_atoms_id(i_atom,i_reac))
                  mass_species = mass_species+atom_mass*DBLE(reac_stoich(i_atom,i_reac))
               END IF
            END DO
            massfracs_reactants(id_reactants(i_reac,2)) = n_spec(i_reac) * mass_species
            if (i_reac <= N_gas) then
               MMW = MMW + massfracs_reactants(id_reactants(i_reac,2))/mass_species
               msum = msum + massfracs_reactants(id_reactants(i_reac,2))
            end if
         END DO
         massfracs_reactants = massfracs_reactants / SUM(massfracs_reactants)

         ! Mean molecular weight is only calculated as mean molecular gas weight per all gas species
         MMW = MMW/msum
         MMW = 1d0/MMW

         rho = (press*1d6)/kB/temp*MMW
         MMW = MMW/amu

      END IF

   end subroutine ec_COMP_EQU_CHEM

   !> Initialize all abundances with uniform abundances for gas and 0 for condensates
   subroutine ec_INIT_ALL_VALS(N_atoms_use,N_reac,n,n_spec,pi_atom)
      !! I/O:
      integer, intent(in)           :: N_atoms_use, N_reac
      DOUBLE PRECISION, intent(out) :: n ! Moles of gas particles per total mass of mixture in kg
      DOUBLE PRECISION, intent(out) :: n_spec(N_reac) ! Moles of species per total mass of mixture in kg
      DOUBLE PRECISION, intent(out) :: pi_atom(N_atoms_use) ! Lagrangian multipliers for the atomic species divided
      ! by (R*T)

      !! Internal:
      INTEGER                      :: i_reac

      n = 0.1d0
      n_spec = 0d0
      pi_atom = 0d0
      DO i_reac = 1, N_gas
         n_spec(i_reac) = n/DBLE(N_gas)
         IF (remove_ions) THEN
            IF(reac_ion(i_reac)) THEN
               n_spec(i_reac) = 0d0
            END IF
         END IF
      END DO

   end subroutine ec_INIT_ALL_VALS

   !> Selects which solid to include next
   subroutine ec_INCLUDE_WHICH_SOLID(N_atoms_use,N_reac,pi_atom,H_0,S_0,temp, &
      n_spec,solid_inclu,neg_cond,dgdnj,remove_cond,inc_next)
      !! I/O
      integer, intent(in)           :: N_atoms_use, N_reac
      DOUBLE PRECISION, intent(in)  :: pi_atom(N_atoms_use)
      DOUBLE PRECISION, intent(in)  :: H_0(N_reac), S_0(N_reac), temp
      LOGICAL, intent(in)           :: solid_inclu(N_reac-N_gas), neg_cond(N_reac-N_gas)

      double precision, intent(inout)  :: n_spec(N_reac), dgdnj(N_reac-N_gas)
      logical, intent(out)          :: remove_cond
      integer, intent(out)          :: inc_next

      !! Internal:
      DOUBLE PRECISION             :: a(N_reactants,N_atoms_use)
      DOUBLE PRECISION             :: mu(N_reactants), minval_inc
      INTEGER                      :: i_atom, i_reac, i_ratom, remove_count, remove_id

      !f2py integer, intent(aux) :: N_gas

      remove_count = 0
      remove_id = 0
      remove_cond = .false.
      DO i_reac = N_gas+1, N_reactants
         IF (n_spec(i_reac) < 0d0) THEN
            ! n_spec(i_reac) = 0d0
            remove_count = remove_count + 1
            remove_id = i_reac
            ! remove_cond = .TRUE.
         END IF
      END DO

      if (remove_count >= 2) then
         inc_next = remove_id
         remove_cond = .true.
      else

         inc_next = -1
         minval_inc = 2d0

         ! Set up a_ij
         a = 0d0
         DO i_atom = 1, N_atoms_use
            DO i_reac = 1, N_reactants
               IF (remove_ions) THEN
                  IF (reac_ion(i_reac)) THEN
                     CYCLE
                  END IF
               END IF
               DO i_ratom = 1, 5
                  IF (reac_atoms_id(i_ratom, i_reac)>0 .and. id_atoms(i_atom) == reac_atoms_id(i_ratom, i_reac)) then
                     a(i_reac,i_atom) = reac_stoich(i_ratom,i_reac)*mol
                  END IF
               END DO
            END DO
         END DO

         ! EVAL Eq. 3.7 in McBride Manual
         DO i_reac = N_gas+1, N_reactants
            mu(i_reac) = H_0(i_reac) - temp*S_0(i_reac)
            dgdnj(i_reac-N_gas) = mu(i_reac)/R/temp - SUM(a(i_reac,1:N_atoms_use)*pi_atom)
         END DO

         DO i_reac = N_gas+1, N_reactants
            IF ((dgdnj(i_reac-N_gas) < 0d0) .AND. (.NOT. solid_inclu(i_reac-N_gas))) THEN
               IF (((dgdnj(i_reac-N_gas) < minval_inc) .AND. (.NOT. neg_cond(i_reac-N_gas)) .AND. &
               temp <= thermo_data_temps(2,thermo_data_n_intervs(i_reac),i_reac)) .AND. &
               (temp >= thermo_data_temps(1,1,i_reac))) THEN
                  minval_inc = dgdnj(i_reac-N_gas)
                  inc_next = i_reac
               END IF
            END IF
         END DO

         if (inc_next==-1 .and. remove_count==1) then
            inc_next = remove_id
            remove_cond = .true.
         end if

      end if

   end subroutine ec_INCLUDE_WHICH_SOLID

   !> Initialize one condensate abundance, as if the most molecules condensed
   subroutine ec_INIT_COND_VALS(N_atoms_use, N_reac, molfracs_atoms, i_cond, n_spec)

      integer, intent(in)              :: N_atoms_use, i_cond, N_reac
      double precision, intent(in)     :: molfracs_atoms(N_atoms_use)
      double precision, intent(inout)  :: n_spec(N_reac)

      integer           :: i_ratom, i_atom
      double precision  :: min_molfrac, stoich_molfrac

      min_molfrac = -1
      do i_ratom = 1, 5
         if (reac_atoms_id(i_ratom,i_cond) > 0) then
            do i_atom = 1, N_atoms_use
               if (id_atoms(i_atom) == reac_atoms_id(i_ratom,i_cond)) then
                  stoich_molfrac = molfracs_atoms(i_atom) / reac_stoich(i_ratom,i_cond)
                  if (min_molfrac==-1 .or. stoich_molfrac < min_molfrac) then
                     min_molfrac = stoich_molfrac
                  end if
               end if
            end do
         end if
      end do

      if (min_molfrac >= 0) then
         n_spec(i_cond) = min_molfrac
      else
         print *, 'ERROR: no data found for the atoms of the given condensate'
      end if

   end subroutine ec_INIT_COND_VALS

   !> Build the small matrix
   subroutine ec_PREP_MATRIX_SHORT(N_atoms_use, N_reac, molfracs_atoms, N_species, press, temp, &
   H_0, S_0, n, n_spec, matrix, vector, solid_indices, N_solids, mu_gas, a_gas)
      !! I/O:
      INTEGER, intent(in)          :: N_atoms_use, N_reac, N_species, N_solids
      INTEGER, intent(in)          :: solid_indices(N_solids)
      DOUBLE PRECISION, intent(in) :: molfracs_atoms(N_atoms_use), press, temp
      DOUBLE PRECISION, intent(in) ::  H_0(N_reac), S_0(N_reac)
      DOUBLE PRECISION, intent(in) :: n ! Moles of gas particles per total mass of mixture in kg
      DOUBLE PRECISION, intent(inout)  :: n_spec(N_reac) ! Moles of species per total mass of mixture in kg

      DOUBLE PRECISION, intent(out):: matrix(N_atoms_use+1+(N_species-N_gas),N_atoms_use+1+(N_species-N_gas))
      ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for
      ! condensed species, the pis and the delta log(n)
      DOUBLE PRECISION, intent(out):: vector(N_atoms_use+1+(N_species-N_gas))
      DOUBLE PRECISION, intent(out):: mu_gas(N_gas), a_gas(N_gas,N_atoms_use)

      !! Internal:
      DOUBLE PRECISION             :: b_0(N_atoms_use), b_0_norm, b(N_atoms_use)
      DOUBLE PRECISION             :: a(N_species,N_atoms_use), mu(N_species)
      INTEGER                      :: i_atom, i_reac, i_ratom, i_atom2

      !f2py integer, intent(aux) :: N_gas

      ! print *, "START ec_PREP_MATRIX_SHORT"

      ! Set up b0
      ! b_0_norm = 0d0
      ! DO i_atom = 1, N_atoms_use
      !    ! call ec_ATOM_MASS(names_atoms(i_atom),mass_atom)
      !    mass_atom = masses_atoms_save(id_atoms(i_atom))
      !    b_0_norm = b_0_norm + mass_atom*molfracs_atoms(i_atom)
      ! END DO
      ! b_0 = molfracs_atoms/b_0_norm
      call ec_b_0(N_atoms_use, molfracs_atoms, b_0_norm, b_0)

      ! Set up a_ij
      a = 0d0
      DO i_atom = 1, N_atoms_use
         ! call uppercase(names_atoms(i_atom),upper_atom_name)
         DO i_reac = 1, N_gas
            IF (remove_ions) THEN
               IF (reac_ion(i_reac)) THEN
                  a(i_reac,1:N_atoms_use) = 0d0
                  CYCLE
               END IF
            END IF
            DO i_ratom = 1, 5
               IF (reac_atoms_id(i_ratom, i_reac)>0 .and. &
               id_atoms(i_atom) == reac_atoms_id(i_ratom, i_reac)) then
                  a(i_reac,i_atom) = reac_stoich(i_ratom,i_reac)*mol
               END IF
            END DO
         END DO
         DO i_reac = N_gas+1, N_species
            DO i_ratom = 1, 5
               IF (reac_atoms_id(i_ratom, solid_indices(i_reac - N_gas))>0 .and. &
               id_atoms(i_atom) == reac_atoms_id(i_ratom, solid_indices(i_reac - N_gas))) then
                  a(i_reac,i_atom) = reac_stoich(i_ratom,solid_indices(i_reac-N_gas))*mol
                  ! print *, i_ratom, i_reac, i_atom
               END IF
            END DO
         END DO
      END DO

      ! Set up mu_j
      DO i_reac = 1, N_species
         IF (remove_ions) THEN
            IF (reac_ion(i_reac)) THEN
               mu(i_reac) = 0d0
               CYCLE
            END IF
         END IF
         ! Taken from Venot et al. (2012), in comparison with McBride 1996.
         IF (i_reac <= N_gas) THEN
            mu(i_reac) = H_0(i_reac) - temp*S_0(i_reac)

            IF (n_spec(i_reac) > 1d-290) THEN
               mu(i_reac) = mu(i_reac) + R*temp*log(n_spec(i_reac)/n)+R*temp*log(press)
            ELSE
               IF (verbose) THEN
                  write(*,*) 'n_spec(i_reac) == 0 for '//trim(adjustl(names_reactants(i_reac)))// &
                  ' set to 1d-13 and try again.'
               END IF
               call RANDOM_NUMBER(n_spec(i_reac))
               n_spec(i_reac) = n_spec(i_reac)*1d-13
               mu(i_reac) = mu(i_reac) + R*temp*log(n_spec(i_reac)/n)+R*temp*log(press)
            END IF
         ELSE
            mu(i_reac) = H_0(solid_indices(i_reac-N_gas)) - temp*S_0(solid_indices(i_reac-N_gas))
         END IF
      END DO

      a_gas = a(1:N_gas,1:N_atoms_use)
      mu_gas = mu(1:N_gas)

      ! MATRIX SETUP
      matrix = 0d0

      ! Set up the matrix for the N_atoms equations (Eq. 2.24)
      DO i_atom = 1, N_atoms_use
         DO i_atom2 = 1, N_atoms_use
            DO i_reac = 1, N_gas
               ! IF (remove_ions) THEN
               !    IF (reac_ion(i_reac)) THEN
               !       CYCLE
               !    END IF
               ! END IF
               if (.not. remove_ions .or. .not. reac_ion(i_reac)) then
                  matrix(i_atom,i_atom2) = matrix(i_atom,i_atom2) + &
                  a(i_reac,i_atom)*a(i_reac,i_atom2)*n_spec(i_reac)
               end if
            END DO
         END DO

         DO i_reac = 1, N_gas
            ! IF (remove_ions) THEN
            !    IF (reac_ion(i_reac)) THEN
            !       CYCLE
            !    END IF
            ! END IF
            if (.not. remove_ions .or. .not. reac_ion(i_reac)) then
               matrix(i_atom,N_atoms_use+1) = matrix(i_atom,N_atoms_use+1) + &
               a(i_reac,i_atom)*n_spec(i_reac)
            end if
         END DO

         IF (N_gas < N_species) THEN
            DO i_reac = N_gas+1, N_species
               matrix(i_atom,N_atoms_use+1+i_reac-N_gas) = a(i_reac,i_atom)
            END DO
         END IF

      END DO

      ! Set up the matrix for the equation (Eq. 2.26)
      DO i_atom = 1, N_atoms_use
         DO i_reac = 1, N_gas
            IF (remove_ions) THEN
               IF (reac_ion(i_reac)) THEN
                  CYCLE
               END IF
            END IF
            matrix(N_atoms_use+1,i_atom) = matrix(N_atoms_use+1,i_atom) + &
            a(i_reac,i_atom)*n_spec(i_reac)
         END DO
      END DO

      DO i_reac = 1, N_gas
         IF (remove_ions) THEN
            IF (reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         matrix(N_atoms_use+1,N_atoms_use+1) = matrix(N_atoms_use+1,N_atoms_use+1) + n_spec(i_reac) !!
      END DO
      matrix(N_atoms_use+1,N_atoms_use+1) = matrix(N_atoms_use+1,N_atoms_use+1) - n

      ! Set up the matrix for the (N_reactants-N_gas) equations (Eq. 2.25)

      IF (N_gas < N_species) THEN
         DO i_reac = N_gas+1, N_species
            DO i_atom = 1, N_atoms_use
               matrix(N_atoms_use+1+i_reac-N_gas,i_atom) = a(i_reac,i_atom)
            END DO
         END DO
      END IF

      ! VECTOR SETUP
      !vector(N_atoms+1+(N_reactants-N_gas))
      vector = 0d0

      ! (Eq. 2.25)
      IF (N_gas < N_species) THEN
         vector(N_atoms_use+2:N_atoms_use+1+(N_species-N_gas)) = mu(N_gas+1:N_species)/R/temp
      END IF

      ! (Eq. 2.24)
      b = 0d0
      DO i_atom = 1, N_atoms_use
         DO i_reac = 1, N_gas
            IF (remove_ions) THEN
               IF (reac_ion(i_reac)) THEN
                  CYCLE
               END IF
            END IF
            b(i_atom) = b(i_atom) + a(i_reac,i_atom)*n_spec(i_reac)
         END DO
         DO i_reac = N_gas+1, N_species
            b(i_atom) = b(i_atom) + a(i_reac,i_atom)*n_spec(solid_indices(i_reac-N_gas))
         END DO
      END DO
      vector(1:N_atoms_use) = b_0 - b
      DO i_reac = 1, N_gas
         IF (remove_ions) THEN
            IF (reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         vector(1:N_atoms_use) = vector(1:N_atoms_use) + &
         a(i_reac,1:N_atoms_use)*n_spec(i_reac)*mu(i_reac)/R/temp
      END DO

      ! (Eq. 2.26)
      vector(N_atoms_use+1) = n - SUM(n_spec(1:N_gas)) + SUM(n_spec(1:N_gas)*mu(1:N_gas))/R/temp

   end subroutine ec_PREP_MATRIX_SHORT



   !> Return the result of one step of computation with small matrix(problem: AX=B)
   subroutine ec_UPDATE_ABUNDS_SHORT(N_atoms_use,N_reac,N_species,solution_vector,n_spec,pi_atom,&
   n,converged,solid_indices,N_solids,mu_gas,a_gas,temp,molfracs_atoms,n_spec_old)

      !! I/O:
      INTEGER, intent(in)          :: N_atoms_use, N_reac, N_species, N_solids
      INTEGER, intent(in)          :: solid_indices(N_solids)
      DOUBLE PRECISION, intent(in) :: solution_vector(N_atoms_use+1+(N_species-N_gas))
      DOUBLE PRECISION, intent(inout)  :: n ! Moles of gas particles per total mass of mixture in kg
      DOUBLE PRECISION, intent(inout)  :: n_spec(N_reac) ! Moles of species per total mass of mixture in kg
      DOUBLE PRECISION, intent(in) :: n_spec_old(N_reac) ! Moles of species per total mass of mixture in kg
      DOUBLE PRECISION, intent(inout)  :: pi_atom(N_atoms_use) ! Lagrangian multipliers for the atomic species divided
      ! by (R*T)
      LOGICAL, intent(out)          :: converged
      DOUBLE PRECISION, intent(in) :: mu_gas(N_gas), a_gas(N_gas,N_atoms_use), temp

      !! Internal:
      INTEGER                      :: i_reac
      INTEGER, save                :: n_done = 0
      DOUBLE PRECISION             :: lambda, lambda1, lambda2
      DOUBLE PRECISION, parameter  :: SIZE = 18.420681
      LOGICAL                      :: gas_good, solids_good, total_good
      DOUBLE PRECISION             :: delta_n_gas(N_gas)

      ! IONS
      INTEGER                      :: i_ion, i_stoich
      DOUBLE PRECISION             :: pi_ion, pi_ion_norm

      ! MASS BALANCE CHECKS
      DOUBLE PRECISION             :: b_0(N_atoms_use), b_0_norm, pi_atom_old(N_atoms_use)
      DOUBLE PRECISION             :: a(N_species,N_atoms_use), mval_mass_good
      INTEGER                      :: i_atom, i_ratom
      LOGICAL                      :: mass_good, pi_good
      DOUBLE PRECISION             :: molfracs_atoms(N_atoms_use)
      DOUBLE PRECISION             :: change

      !f2py integer, intent(aux) :: N_gas

      ! print *, "START ec_UPDATE_ABUNDS_SHORT"

      ! Get delta_n_gas, following Eq. 2.18:
      DO i_reac = 1, N_gas
         IF (remove_ions) THEN
            IF (reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         delta_n_gas(i_reac) = SUM(a_gas(i_reac,1:N_atoms_use)*solution_vector(1:N_atoms_use)) + &
         solution_vector(N_atoms_use+1) - mu_gas(i_reac)/R/temp
      END DO

      ! Calculate correction factors as described in Section 3.3 of the McBride Manual
      lambda1 = 9d99
      lambda2 = 9d99
      DO i_reac = 1, N_gas
         IF (remove_ions) THEN
            IF (reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         IF (LOG(n_spec(i_reac)/n) > -SIZE) THEN
            lambda1 = MIN(lambda1,2d0/(MAX(5d0*ABS(solution_vector(N_atoms_use+1)), &
            ABS(delta_n_gas(i_reac)))))
         ELSE IF ((LOG(n_spec(i_reac)/n) <= -SIZE) .AND. (delta_n_gas(i_reac) >= 0d0)) THEN
            lambda2 = MIN(lambda2,ABS((-LOG(n_spec(i_reac)/n)-9.2103404)/ &
            (delta_n_gas(i_reac)-solution_vector(N_atoms_use+1))))
         END IF
      END DO
      lambda = MIN(1d0,lambda1,lambda2)

      DO i_reac = 1, N_gas
         IF (remove_ions) THEN
            IF (reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         n_spec(i_reac) = n_spec(i_reac)*exp(lambda*delta_n_gas(i_reac))
      END DO

      IF (N_gas < N_species) THEN
         DO i_reac = N_gas+1, N_species
            change = lambda*solution_vector(N_atoms_use+1+i_reac-N_gas)
            ! if (2*abs(change) < n_spec(solid_indices(i_reac-N_gas)) .OR. n_spec(solid_indices(i_reac-N_gas)) < tiny(0d0)) then
            !    n_spec(solid_indices(i_reac-N_gas)) = n_spec(solid_indices(i_reac-N_gas)) + change
            ! else
            !    n_spec(solid_indices(i_reac-N_gas)) = (1d0 + sign(0.5d0,change)) * n_spec(solid_indices(i_reac-N_gas))
            ! end if
            n_spec(solid_indices(i_reac-N_gas)) = n_spec(solid_indices(i_reac-N_gas)) + change
         END DO
      END IF
      pi_atom_old = pi_atom
      pi_atom = solution_vector(1:N_atoms_use)
      n = n*exp(lambda*solution_vector(N_atoms_use+1))

      gas_good = .TRUE.
      DO i_reac = 1, N_gas
         IF (remove_ions) THEN
            IF (reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         IF (n_spec(i_reac)*ABS(delta_n_gas(i_reac))/SUM(n_spec) > 0.5d-5) THEN
            gas_good = .FALSE.
         END IF
      END DO
      solids_good = .TRUE.
      IF (N_gas < N_species) THEN
         DO i_reac = N_gas+1, N_species
            IF (ABS(solution_vector(N_atoms_use+1+i_reac-N_gas))/SUM(n_spec) > 0.5d-5) THEN
               solids_good = .FALSE.
            END IF
         END DO
      END IF
      total_good = .TRUE.
      IF (n*ABS(solution_vector(N_atoms_use+1))/SUM(n_spec) > 05.d-5) THEN
         total_good = .FALSE.
      END IF

      !!!-----------------------

      mass_good = .TRUE.
      pi_good = .TRUE.

      ! Set up b0
      ! b_0_norm = 0d0
      ! DO i_atom = 1, N_atoms_use
      !    ! call ec_ATOM_MASS(names_atoms(i_atom),mass_atom)
      !    mass_atom = masses_atoms_save(id_atoms(i_atom))
      !    b_0_norm = b_0_norm + mass_atom*molfracs_atoms(i_atom)
      ! END DO
      ! b_0 = molfracs_atoms/b_0_norm
      call ec_b_0(N_atoms_use, molfracs_atoms, b_0_norm, b_0)

      ! Set up a_ij
      a = 0d0
      DO i_atom = 1, N_atoms_use
         DO i_reac = 1, N_gas
            IF (remove_ions) THEN
               IF (reac_ion(i_reac)) THEN
                  a(i_reac,1:N_atoms_use) = 0d0
                  CYCLE
               END IF
            END IF
            DO i_ratom = 1, 5
               IF (reac_atoms_id(i_ratom, i_reac)>0 .and. id_atoms(i_atom) == reac_atoms_id(i_ratom, i_reac)) then
                  a(i_reac,i_atom) = reac_stoich(i_ratom,i_reac)*mol
               END IF
            END DO
         END DO
         DO i_reac = N_gas+1, N_species
            DO i_ratom = 1, 5
               IF (reac_atoms_id(i_ratom, solid_indices(i_reac - N_gas))>0 .and. &
               id_atoms(i_atom) == reac_atoms_id(i_ratom, solid_indices(i_reac - N_gas))) then
                  a(i_reac,i_atom) = reac_stoich(i_ratom,solid_indices(i_reac-N_gas))*mol
               END IF
            END DO
         END DO
      END DO

      mval_mass_good = MAXVAL(b_0)*1d-2
      DO i_atom = 1, N_atoms_use
         IF ((abs(b_0(i_atom)-sum(a(1:N_species,i_atom)*n_spec(1:N_species))) > mval_mass_good) .AND. (b_0(i_atom) > 1d-6)) THEN
            mass_good = .FALSE.
         END IF
      END DO

      DO i_atom = 1, N_atoms_use
         IF (abs((pi_atom_old(i_atom)-pi_atom(i_atom))/pi_atom(i_atom)) > 1d-3) THEN
            pi_good = .FALSE.
         END IF
      END DO

      IF ((.NOT. mass_good) .OR. (.NOT. pi_good)) THEN
         mass_good = .TRUE.
         pi_good = .TRUE.
         DO i_reac = 1, N_reactants
            IF (ABS(n_spec(i_reac)-n_spec_old(i_reac)) > 1d-10) THEN
               mass_good = .FALSE.
               pi_good = .FALSE.
            END IF
         END DO
      END IF

      !!!-------------------

      ! ION CONVERGENCE?

      IF (ions .AND. (.NOT. remove_ions)) THEN

         ! DO THE MAGIC THEY DO IN SECT. 3.7 in McBride

         pi_ion = 0d0
         pi_ion_norm = 0d0
         DO i_reac = 1, N_species
            DO i_stoich = 1, 5
               ! IF (trim(adjustl(reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
               IF (reac_atoms_id(i_stoich,i_reac) == 1) THEN
                  pi_ion = pi_ion - n_spec(i_reac)*reac_stoich(i_stoich,i_reac)
                  pi_ion_norm = pi_ion_norm + n_spec(i_reac)*reac_stoich(i_stoich,i_reac)**2d0
                  EXIT
               END IF
            END DO
         END DO

         pi_ion = pi_ion / pi_ion_norm

         IF (ABS(pi_ion) > 1d-4) THEN
            DO i_ion = 1, 80
               DO i_reac = 1, N_species
                  DO i_stoich = 1, 5
                     ! IF (trim(adjustl(reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
                     IF (reac_atoms_id(i_stoich,i_reac) == 1) THEN
                        n_spec(i_reac) = n_spec(i_reac)*exp(reac_stoich(i_stoich,i_reac)*pi_ion)
                        EXIT
                     END IF
                  END DO
               END DO

               pi_ion = 0d0
               pi_ion_norm = 0d0
               DO i_reac = 1, N_species
                  DO i_stoich = 1, 5
                     ! IF (trim(adjustl(reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
                     IF (reac_atoms_id(i_stoich,i_reac) == 1) THEN
                        pi_ion = pi_ion - n_spec(i_reac)*reac_stoich(i_stoich,i_reac)
                        pi_ion_norm = pi_ion_norm + n_spec(i_reac)*reac_stoich(i_stoich,i_reac)**2d0
                        EXIT
                     END IF
                  END DO
               END DO

            END DO
         END IF

         IF (((gas_good .AND. solids_good) .AND. total_good) .AND. (ABS(pi_ion) <= 1d-4)) THEN
            IF (mass_good .AND. pi_good) THEN
               converged = .TRUE.
            END IF
         END IF

      ELSE

         IF ((gas_good .AND. solids_good) .AND. total_good) THEN
            IF (mass_good .AND. pi_good) THEN
               converged = .TRUE.
            END IF
         END IF

      END IF

      n_done = n_done + 1

   end subroutine ec_UPDATE_ABUNDS_SHORT

   !> Build the big matrix
   subroutine ec_PREP_MATRIX_LONG(N_atoms_use, N_reac, molfracs_atoms, N_species,press,temp, &
   H_0, S_0, n, n_spec, matrix, vector, solid_indices, N_solids)

      !! I/O:
      INTEGER, intent(in)          :: N_atoms_use, N_reac, N_species, N_solids
      INTEGER, intent(in)          :: solid_indices(N_solids)
      DOUBLE PRECISION, intent(in) :: molfracs_atoms(N_atoms_use), press, temp
      DOUBLE PRECISION, intent(in) :: H_0(N_reac), S_0(N_reac)
      DOUBLE PRECISION, intent(in) :: n ! Moles of gas particles per total mass of mixture in kg
      DOUBLE PRECISION, intent(inout)  :: n_spec(N_reac) ! Moles of species per total mass of mixture in kg
      ! DOUBLE PRECISION, intent(in) :: pi_atom(N_atoms) ! Lagrangian multipliers for the atomic species divided
      ! ! by (R*T)
      DOUBLE PRECISION, intent(out):: matrix(N_species+N_atoms_use+1,N_species+N_atoms_use+1)
      ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for
      ! condensed species, the pis and the delta log(n)
      DOUBLE PRECISION, intent(out):: vector(N_species+N_atoms_use+1)

      !! Internal:
      DOUBLE PRECISION             :: b_0(N_atoms_use), b_0_norm, b(N_atoms_use)
      DOUBLE PRECISION             :: a(N_species,N_atoms_use), mu(N_species)
      INTEGER                      :: i_atom, i_reac, i_ratom

      ! Set up b0
      ! b_0_norm = 0d0
      ! DO i_atom = 1, N_atoms_use
      !    ! call ec_ATOM_MASS(names_atoms(i_atom),mass_atom)
      !    if (id_atoms(i_atom) > 0) then
      !       mass_atom = masses_atoms_save(id_atoms(i_atom))
      !       b_0_norm = b_0_norm + mass_atom*molfracs_atoms(i_atom)
      !    end if
      ! END DO
      ! b_0 = molfracs_atoms/b_0_norm
      call ec_b_0(N_atoms_use, molfracs_atoms, b_0_norm, b_0)

      ! Set up a_ij
      a = 0d0
      DO i_atom = 1, N_atoms_use
         ! call uppercase(names_atoms(i_atom),upper_atom_name)
         DO i_reac = 1, N_gas
            IF (remove_ions) THEN
               IF (reac_ion(i_reac)) THEN
                  a(i_reac,1:N_atoms_use) = 0d0
                  CYCLE
               END IF
            END IF
            DO i_ratom = 1, 5
               ! call uppercase(reac_atoms_names(i_ratom,i_reac),upper_ratom_name)
               ! IF (trim(adjustl(upper_atom_name)) .EQ. trim(adjustl(upper_ratom_name))) THEN
               IF (reac_atoms_id(i_ratom,i_reac)>0 .and. id_atoms(i_atom) == reac_atoms_id(i_ratom, i_reac)) then
                  a(i_reac,i_atom) = reac_stoich(i_ratom,i_reac)*mol
               END IF
            END DO
         END DO
         DO i_reac = N_gas+1, N_species
            DO i_ratom = 1, 5
               ! call uppercase(reac_atoms_names(i_ratom,solid_indices(i_reac-N_gas)),upper_ratom_name)
               ! IF (trim(adjustl(upper_atom_name)) .EQ. trim(adjustl(upper_ratom_name))) THEN
               IF (reac_atoms_id(i_ratom, solid_indices(i_reac - N_gas))>0 .and.&
               id_atoms(i_atom) == reac_atoms_id(i_ratom, solid_indices(i_reac - N_gas))) then
                  a(i_reac,i_atom) = reac_stoich(i_ratom,solid_indices(i_reac-N_gas))*mol
               END IF
            END DO
         END DO
      END DO

      ! Set up mu_j
      DO i_reac = 1, N_species
         IF (remove_ions) THEN
            IF (reac_ion(i_reac)) THEN
               mu(i_reac) = 0d0
               CYCLE
            END IF
         END IF
         ! Taken from Venot et al. (2012), in comparison with McBride 1996.
         IF (i_reac <= N_gas) THEN
            mu(i_reac) = H_0(i_reac) - temp*S_0(i_reac)

            IF (n_spec(i_reac) > 1d-290) THEN
               mu(i_reac) = mu(i_reac) + R*temp*log(n_spec(i_reac)/n)+R*temp*log(press)
            ELSE
               IF (verbose) THEN
                  write(*,*) 'n_spec(i_reac) == 0 for '//trim(adjustl(names_reactants(i_reac)))// &
                  ' set to 1d-13 and try again.'
               END IF
               call RANDOM_NUMBER(n_spec(i_reac))
               n_spec(i_reac) = n_spec(i_reac)*1d-13
               mu(i_reac) = mu(i_reac) + R*temp*log(n_spec(i_reac)/n)+R*temp*log(press)
            END IF

         ELSE
            mu(i_reac) = H_0(solid_indices(i_reac-N_gas)) - temp*S_0(solid_indices(i_reac-N_gas))
         END IF
      END DO

      ! MATRIX SETUP
      matrix = 0d0
      ! Set up the matrix for the N_gas equations (Eq. 2.18)
      DO i_reac = 1, N_gas
         IF (remove_ions) THEN
            IF (reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         matrix(i_reac,i_reac) = 1d0
         DO i_atom = 1, N_atoms_use
            matrix(i_reac,N_species+i_atom) = -a(i_reac,i_atom)
         END DO
         matrix(i_reac,N_species+N_atoms_use+1) = -1d0
      END DO

      ! Set up the matrix for the N_reactants-N_gas equations (Eq. 2.19)
      IF (N_gas < N_species) THEN
         DO i_reac = N_gas+1, N_species
            DO i_atom = 1, N_atoms_use
               matrix(i_reac,N_species+i_atom) = -a(i_reac,i_atom)
            END DO
         END DO
      END IF

      ! Set up the matrix for the N_atom equations (Eq. 2.20)
      DO i_atom = 1, N_atoms_use
         DO i_reac = 1, N_gas
            IF (remove_ions) THEN
               IF (reac_ion(i_reac)) THEN
                  CYCLE
               END IF
            END IF
            matrix(N_species+i_atom,i_reac) = a(i_reac,i_atom)*n_spec(i_reac)
         END DO
         IF (N_gas < N_species) THEN
            DO i_reac = N_gas+1, N_species
               matrix(N_species+i_atom,i_reac) = a(i_reac,i_atom)
            END DO
         END IF
      END DO

      ! Set up the matrix for the last equation (Eq. 2.21)
      DO i_reac = 1, N_gas
         IF (remove_ions) THEN
            IF (reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         matrix(N_species+N_atoms_use+1,i_reac) = n_spec(i_reac)
      END DO
      matrix(N_species+N_atoms_use+1,N_species+N_atoms_use+1) = -n

      ! VECTOR SETUP
      !vector(N_reactants+N_atoms+1)
      vector = 0d0

      DO i_reac = 1, N_gas
         IF (remove_ions) THEN
            IF (reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         vector(i_reac) = -mu(i_reac)/R/temp ! (Eq. 2.18)
      END DO

      IF (N_gas < N_species) THEN
         vector(N_gas+1:N_species) = -mu(N_gas+1:N_species)/R/temp ! (Eq. 2.19)
      END IF

      b = 0d0
      DO i_atom = 1, N_atoms_use
         DO i_reac = 1, N_gas
            IF (remove_ions) THEN
               IF (reac_ion(i_reac)) THEN
                  CYCLE
               END IF
            END IF
            b(i_atom) = b(i_atom) + a(i_reac,i_atom)*n_spec(i_reac)
         END DO
         DO i_reac = N_gas+1, N_species
            b(i_atom) = b(i_atom) + a(i_reac,i_atom)*n_spec(solid_indices(i_reac-N_gas))
         END DO
      END DO
      vector(N_species+1:N_species+N_atoms_use) = b_0 - b ! (Eq. 2.20)

      vector(N_species+N_atoms_use+1) = n - SUM(n_spec(1:N_gas)) ! (Eq. 2.21)

   end subroutine ec_PREP_MATRIX_LONG

   !> Return the result of one step of computation with big matrix(problem: AX=B)
   subroutine ec_UPDATE_ABUNDS_LONG(N_atoms_use,N_reac,N_species,solution_vector,n_spec,pi_atom,&
   n,converged,solid_indices,N_solids,molfracs_atoms,n_spec_old)

      !! I/O:
      INTEGER, intent(in)          :: N_atoms_use, N_reac, N_species, N_solids
      INTEGER , intent(in)         :: solid_indices(N_solids)
      DOUBLE PRECISION, intent(in) :: solution_vector(N_species+N_atoms_use+1)
      DOUBLE PRECISION, intent(inout)  :: n ! Moles of gas particles per total mass of mixture in kg
      DOUBLE PRECISION, intent(inout)  :: n_spec(N_reac) ! Moles of species per total mass of mixture in kg
      DOUBLE PRECISION, intent(in)     :: n_spec_old(N_reac) ! Moles of species per total mass of mixture in kg
      DOUBLE PRECISION, intent(inout)  :: pi_atom(N_atoms_use) ! Lagrangian multipliers for the atomic species divided
      ! by (R*T)
      LOGICAL , intent(out)         :: converged

      !! Internal:
      INTEGER                      :: i_reac
      INTEGER, save                :: n_done = 0
      DOUBLE PRECISION             :: lambda, lambda1, lambda2
      DOUBLE PRECISION, parameter  :: SIZE = 18.420681
      LOGICAL                      :: gas_good, solids_good, total_good

      ! IONS
      INTEGER                      :: i_ion, i_stoich
      DOUBLE PRECISION             :: pi_ion, pi_ion_norm

      ! MASS BALANCE CHECKS
      DOUBLE PRECISION             :: b_0(N_atoms_use), b_0_norm, pi_atom_old(N_atoms_use)
      DOUBLE PRECISION             :: a(N_species,N_atoms_use), mval_mass_good
      INTEGER                      :: i_atom, i_ratom
      LOGICAL                      :: mass_good, pi_good
      DOUBLE PRECISION             :: molfracs_atoms(N_atoms_use)
      DOUBLE precision             :: change


      ! Calculate correction factors as described in Section 3.3 of the McBride Manual
      lambda1 = 9d99
      lambda2 = 9d99
      DO i_reac = 1, N_gas
         IF (remove_ions) THEN
            IF (reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         IF (LOG(n_spec(i_reac)/n) > -SIZE) THEN
            lambda1 = MIN(lambda1,2d0/(MAX(5d0*ABS(solution_vector(N_species+N_atoms_use+1)), &
            ABS(solution_vector(i_reac)))))
         ELSE IF ((LOG(n_spec(i_reac)/n) <= -SIZE) .AND. (solution_vector(i_reac) >= 0d0)) THEN
            lambda2 = MIN(lambda2,ABS((-LOG(n_spec(i_reac)/n)-9.2103404)/ &
            (solution_vector(i_reac)-solution_vector(N_species+N_atoms_use+1))))
         END IF
      END DO
      lambda = MIN(1d0,lambda1,lambda2)

      DO i_reac = 1, N_gas
         IF (remove_ions) THEN
            IF (reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         n_spec(i_reac) = n_spec(i_reac)*exp(lambda*solution_vector(i_reac))
      END DO

      IF (N_gas < N_species) THEN
         DO i_reac = N_gas+1, N_species
            change = lambda*solution_vector(N_atoms_use+1+i_reac-N_gas)
            ! if (2*abs(change) < n_spec(solid_indices(i_reac-N_gas)) .OR. n_spec(solid_indices(i_reac-N_gas)) < tiny(0d0)) then
            !    n_spec(solid_indices(i_reac-N_gas)) = n_spec(solid_indices(i_reac-N_gas)) + changepi_atom
            ! else
            !    n_spec(solid_indices(i_reac-N_gas)) = (1d0 + sign(0.5d0,change)) * n_spec(solid_indices(i_reac-N_gas))
            ! end if
            n_spec(solid_indices(i_reac-N_gas)) = n_spec(solid_indices(i_reac-N_gas)) + change
         END DO
      END IF
      pi_atom_old = pi_atom
      pi_atom = solution_vector(N_species+1:N_species+N_atoms_use)
      n = n*exp(lambda*solution_vector(N_species+N_atoms_use+1))

      gas_good = .TRUE.
      DO i_reac = 1, N_gas
         IF (remove_ions) THEN
            IF (reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         IF (n_spec(i_reac)*ABS(solution_vector(i_reac))/SUM(n_spec) > 0.5d-5) THEN
            gas_good = .FALSE.
         END IF
      END DO
      solids_good = .TRUE.
      IF (N_gas < N_species) THEN
         DO i_reac = N_gas+1, N_species
            IF (ABS(solution_vector(i_reac))/SUM(n_spec) > 0.5d-5) THEN
               solids_good = .FALSE.
            END IF
         END DO
      END IF
      total_good = .TRUE.
      IF (n*ABS(solution_vector(N_species+N_atoms_use+1))/SUM(n_spec) > 05.d-5) THEN
         total_good = .FALSE.
      END IF

      !!!-----------------------

      mass_good = .TRUE.
      pi_good = .TRUE.


      ! Set up b0
      ! b_0_norm = 0d0
      ! DO i_atom = 1, N_atoms_use
      !    ! call ec_ATOM_MASS(names_atoms(i_atom),mass_atom)
      !    mass_atom = masses_atoms_save(id_atoms(i_atom))
      !    b_0_norm = b_0_norm + mass_atom*molfracs_atoms(i_atom)
      ! END DO
      ! b_0 = molfracs_atoms/b_0_norm
      call ec_b_0(N_atoms_use, molfracs_atoms, b_0_norm, b_0)

      ! Set up a_ij
      a = 0d0
      DO i_atom = 1, N_atoms_use
         ! call uppercase(names_atoms(i_atom),upper_atom_name)
         DO i_reac = 1, N_gas
            IF (remove_ions) THEN
               IF (reac_ion(i_reac)) THEN
                  a(i_reac,1:N_atoms_use) = 0d0
                  CYCLE
               END IF
            END IF
            DO i_ratom = 1, 5
               ! call uppercase(reac_atoms_names(i_ratom,i_reac),upper_ratom_name)
               ! IF (trim(adjustl(upper_atom_name)) .EQ. trim(adjustl(upper_ratom_name))) THEN
               IF (reac_atoms_id(i_ratom, i_reac)>0 .and. id_atoms(i_atom) == reac_atoms_id(i_ratom, i_reac)) then
                  a(i_reac,i_atom) = reac_stoich(i_ratom,i_reac)*mol
               END IF
            END DO
         END DO
         DO i_reac = N_gas+1, N_species
            DO i_ratom = 1, 5
               ! call uppercase(reac_atoms_names(i_ratom,solid_indices(i_reac-N_gas)),upper_ratom_name)
               ! IF (trim(adjustl(upper_atom_name)) .EQ. trim(adjustl(upper_ratom_name))) THEN
               IF (reac_atoms_id(i_ratom, solid_indices(i_reac - N_gas))>0 .and. &
               id_atoms(i_atom) == reac_atoms_id(i_ratom, solid_indices(i_reac - N_gas))) then
                  a(i_reac,i_atom) = reac_stoich(i_ratom,solid_indices(i_reac-N_gas))*mol
               END IF
            END DO
         END DO
      END DO

      mval_mass_good = MAXVAL(b_0)*1d-2
      DO i_atom = 1, N_atoms_use
         IF ((abs(b_0(i_atom)-sum(a(1:N_species,i_atom)*n_spec(1:N_species))) > mval_mass_good) .AND. (b_0(i_atom) > 1d-6)) THEN
            mass_good = .FALSE.
         END IF
      END DO

      DO i_atom = 1, N_atoms_use
         IF (abs((pi_atom_old(i_atom)-pi_atom(i_atom))/pi_atom(i_atom)) > 1d-3) THEN
            pi_good = .FALSE.
         END IF
      END DO

      IF ((.NOT. mass_good) .OR. (.NOT. pi_good)) THEN
         mass_good = .TRUE.
         pi_good = .TRUE.
         DO i_reac = 1, N_reactants
            IF (ABS(n_spec(i_reac)-n_spec_old(i_reac)) > 1d-10) THEN
               mass_good = .FALSE.
               pi_good = .FALSE.
            END IF
         END DO
      END IF

      !!!-------------------

      ! ION CONVERGENCE?

      IF (ions .AND. (.NOT. remove_ions)) THEN

         ! DO THE MAGIC THEY DO IN SECT. 3.7 in McBride

         pi_ion = 0d0
         pi_ion_norm = 0d0
         DO i_reac = 1, N_species
            DO i_stoich = 1, 5
               ! IF (trim(adjustl(reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
               IF (reac_atoms_id(i_stoich,i_reac) == 1) THEN
                  pi_ion = pi_ion - n_spec(i_reac)*reac_stoich(i_stoich,i_reac)
                  pi_ion_norm = pi_ion_norm + n_spec(i_reac)*reac_stoich(i_stoich,i_reac)**2d0
                  EXIT
               END IF
            END DO
         END DO

         pi_ion = pi_ion / pi_ion_norm

         IF (ABS(pi_ion) > 1d-4) THEN
            DO i_ion = 1, 80
               DO i_reac = 1, N_species
                  DO i_stoich = 1, 5
                     ! IF (trim(adjustl(reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
                     IF (reac_atoms_id(i_stoich,i_reac) == 1) THEN
                        n_spec(i_reac) = n_spec(i_reac)*exp(reac_stoich(i_stoich,i_reac)*pi_ion)
                        EXIT
                     END IF
                  END DO
               END DO

               pi_ion = 0d0
               pi_ion_norm = 0d0
               DO i_reac = 1, N_species
                  DO i_stoich = 1, 5
                     ! IF (trim(adjustl(reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
                     IF (reac_atoms_id(i_stoich,i_reac) == 1) THEN
                        pi_ion = pi_ion - n_spec(i_reac)*reac_stoich(i_stoich,i_reac)
                        pi_ion_norm = pi_ion_norm + n_spec(i_reac)*reac_stoich(i_stoich,i_reac)**2d0
                        EXIT
                     END IF
                  END DO
               END DO
            END DO
         END IF

         IF (((gas_good .AND. solids_good) .AND. total_good) .AND. (ABS(pi_ion) <= 1d-4)) THEN
            IF (mass_good .AND. pi_good) THEN
               converged = .TRUE.
            END IF
         END IF

      ELSE

         IF ((gas_good .AND. solids_good) .AND. total_good) THEN
            IF (mass_good .AND. pi_good) THEN
               converged = .TRUE.
            END IF
         END IF

      END IF

      n_done = n_done + 1

   end subroutine ec_UPDATE_ABUNDS_LONG

   !> Computes the adiabatic gradient
   subroutine ec_COMP_ADIABATIC_GRAD(N_atoms_use,N_reac,N_spec_eff,n_spec, &
   n,H_0,C_P_0,solid_indices,N_solids,temp,nabla_ad,gamma2,c_pe)

      !! I/O:
      INTEGER, intent(in)          :: N_atoms_use, N_reac, N_spec_eff, N_solids
      INTEGER, intent(in)          :: solid_indices(N_solids)
      DOUBLE PRECISION, intent(in) :: temp
      DOUBLE PRECISION, intent(in) :: C_P_0(N_reac), H_0(N_reac)
      DOUBLE PRECISION, intent(in) :: n ! Moles of gas particles per total mass of mixture in kg
      DOUBLE PRECISION, intent(in) :: n_spec(N_reac) ! Moles of species per total mass of mixture in kg
      DOUBLE PRECISION, intent(out):: nabla_ad, gamma2, c_pe

      !! Internal:
      DOUBLE PRECISION             :: matrix(N_atoms_use+1+(N_spec_eff-N_gas),N_atoms_use+1+(N_spec_eff-N_gas))
      ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for
      ! condensed species, the pis and the delta log(n)
      DOUBLE PRECISION             :: vector(N_atoms_use+1+(N_spec_eff-N_gas)), &
      solution_vector(N_atoms_use+1+(N_spec_eff-N_gas))
      DOUBLE PRECISION             :: a(N_spec_eff,N_atoms_use)
      INTEGER                      :: i_atom, i_reac, i_ratom, i_atom2

      ! Set up a_ij
      a = 0d0
      DO i_atom = 1, N_atoms_use
         ! call uppercase(names_atoms(i_atom),upper_atom_name)
         DO i_reac = 1, N_gas
            IF (remove_ions) THEN
               IF (reac_ion(i_reac)) THEN
                  a(i_reac,i_atom) = 0d0
                  CYCLE
               END IF
            END IF
            DO i_ratom = 1, 5
               ! call uppercase(reac_atoms_names(i_ratom,i_reac),upper_ratom_name)
               ! IF (trim(adjustl(upper_atom_name)) .EQ. trim(adjustl(upper_ratom_name))) THEN
               IF (reac_atoms_id(i_ratom, i_reac)>0 .and. id_atoms(i_atom) == reac_atoms_id(i_ratom, i_reac)) then
                  a(i_reac,i_atom) = reac_stoich(i_ratom,i_reac)*mol
               END IF
            END DO
         END DO
         DO i_reac = N_gas+1, N_spec_eff
            DO i_ratom = 1, 5
               ! call uppercase(reac_atoms_names(i_ratom,solid_indices(i_reac-N_gas)),upper_ratom_name)
               ! IF (trim(adjustl(upper_atom_name)) .EQ. trim(adjustl(upper_ratom_name))) THEN
               IF (reac_atoms_id(i_ratom, solid_indices(i_reac - N_gas))>0 .and. &
               id_atoms(i_atom) == reac_atoms_id(i_ratom, solid_indices(i_reac - N_gas))) then
                  a(i_reac,i_atom) = reac_stoich(i_ratom,solid_indices(i_reac-N_gas))*mol
               END IF
            END DO
         END DO
      END DO

      matrix = 0d0
      ! Setup matrix, following Eq. 2.56
      DO i_atom = 1, N_atoms_use
         ! First term, LHS
         DO i_atom2 = 1, N_atoms_use
            DO i_reac = 1, N_gas
               IF (remove_ions) THEN
                  IF (reac_ion(i_reac)) THEN
                     CYCLE
                  END IF
               END IF
               matrix(i_atom,i_atom2) = matrix(i_atom,i_atom2) + &
               n_spec(i_reac)*a(i_reac,i_atom2)*a(i_reac,i_atom)
            END DO
         END DO
         ! Second term, LHS
         DO i_reac = N_gas+1, N_spec_eff
            matrix(i_atom,N_atoms_use+1+i_reac-N_gas) = a(i_reac,i_atom)
         END DO
         ! Third term, LHS
         DO i_reac = 1, N_gas
            IF (remove_ions) THEN
               IF (reac_ion(i_reac)) THEN
                  CYCLE
               END IF
            END IF
            matrix(i_atom,N_atoms_use+1) = matrix(i_atom,N_atoms_use+1) + &
            a(i_reac,i_atom)*n_spec(i_reac)
         END DO
      END DO

      ! Setup matrix, following Eq. 2.58
      DO i_atom = 1, N_atoms_use
         DO i_reac = 1, N_gas
            IF (remove_ions) THEN
               IF (reac_ion(i_reac)) THEN
                  CYCLE
               END IF
            END IF
            matrix(N_atoms_use+1,i_atom) = matrix(N_atoms_use+1,i_atom) + &
            a(i_reac,i_atom)*n_spec(i_reac)
         END DO
      END DO

      ! Setup matrix, following Eq. 2.57
      DO i_reac = N_gas+1, N_spec_eff
         DO i_atom = 1, N_atoms_use
            matrix(N_atoms_use+1+i_reac-N_gas,i_atom) = a(i_reac,i_atom)
         END DO
      END DO

      vector = 0d0
      ! Setup the vector, following Eq. 2.56
      DO i_atom = 1, N_atoms_use
         vector(i_atom) = -SUM(a(1:N_gas,i_atom)*n_spec(1:N_gas)*H_0(1:N_gas)) &
         /R/temp
      END DO

      ! Setup the vector, following Eq. 2.58
      vector(N_atoms_use+1) = -SUM(n_spec(1:N_gas)*H_0(1:N_gas))/R/temp

      ! Setup the vector, following Eq. 2.57
      DO i_reac = N_gas+1, N_spec_eff
         vector(N_atoms_use+1+i_reac-N_gas) = -H_0(solid_indices(i_reac-N_gas))/R/temp
      END DO

      ! Solve the system
      call ec_INVERT_MATRIX_SHORT(N_atoms_use+1+N_spec_eff-N_gas,matrix,vector,solution_vector)
      if (error) RETURN

      ! Calculate c_pe, following Eq. 2.59
      c_pe = 0d0
      DO i_atom = 1, N_atoms_use
         c_pe = c_pe + SUM(a(1:N_gas,i_atom)*n_spec(1:N_gas)*H_0(1:N_gas)/R/temp) * &
         solution_vector(i_atom)
      END DO
      DO i_reac = N_gas+1, N_spec_eff
         c_pe = c_pe + H_0(solid_indices(i_reac-N_gas))/R/temp* &
         solution_vector(N_atoms_use+1+i_reac-N_gas) + &
         n_spec(solid_indices(i_reac-N_gas))* &
         C_P_0(solid_indices(i_reac-N_gas))/R
      END DO
      c_pe = c_pe + SUM(n_spec(1:N_gas)*C_P_0(1:N_gas)/R)
      c_pe = c_pe + SUM(n_spec(1:N_gas)*H_0(1:N_gas)/R/temp)* &
      solution_vector(N_atoms_use+1)
      c_pe = c_pe + SUM(n_spec(1:N_gas)*(H_0(1:N_gas)/R/temp)**2d0)
      c_pe = c_pe*R

      ! Calculate nabla_ad, using Eq. 2.50 and Eq. 2.75
      nabla_ad = c_pe/n/R/(1d0+solution_vector(N_atoms_use+1))
      nabla_ad = 1/nabla_ad
      gamma2 = 1d0/(1d0-nabla_ad)

   end subroutine ec_COMP_ADIABATIC_GRAD

   subroutine ec_b_0(N_atoms_use, molfracs_atoms, b_0_norm, b_0)
      integer, intent(in)           :: N_atoms_use
      double precision, intent(in)  :: molfracs_atoms(N_atoms_use)
      double precision, intent(out) :: b_0_norm, b_0(N_atoms_use)

      integer                       :: i_atom
      double precision              :: mass_atom

      ! Set up b0
      b_0_norm = 0d0
      DO i_atom = 1, N_atoms_use
         ! call ec_ATOM_MASS(names_atoms(i_atom),mass_atom)
         mass_atom = masses_atoms_save(id_atoms(i_atom))
         b_0_norm = b_0_norm + mass_atom*molfracs_atoms(i_atom)
      END DO
      b_0 = molfracs_atoms/b_0_norm
   end subroutine ec_b_0

   ! subroutine ec_a(N_atoms_use, a)
   !    integer, intent(in)           :: N_atoms_use
   !    double precision, intent(out) :: a()
   !    a = 0d0
   !    DO i_atom = 1, N_atoms_use
   !       DO i_reac = 1, N_gas
   !          IF (remove_ions) THEN
   !             IF (reac_ion(i_reac)) THEN
   !                a(i_reac,1:N_atoms_use) = 0d0
   !                CYCLE
   !             END IF
   !          END IF
   !          DO i_ratom = 1, 5
   !             IF (reac_atoms_id(i_ratom, i_reac)>0 .and. id_atoms(i_atom) == reac_atoms_id(i_ratom, i_reac)) then
   !                a(i_reac,i_atom) = reac_stoich(i_ratom,i_reac)*mol
   !             END IF
   !          END DO
   !       END DO
   !       DO i_reac = N_gas+1, N_species
   !          DO i_ratom = 1, 5
   !             IF (reac_atoms_id(i_ratom, solid_indices(i_reac - N_gas))>0 .and. &
   !             id_atoms(i_atom) == reac_atoms_id(i_ratom, solid_indices(i_reac - N_gas))) then
   !                a(i_reac,i_atom) = reac_stoich(i_ratom,solid_indices(i_reac-N_gas))*mol
   !             END IF
   !          END DO
   !       END DO
   !    END DO
   ! end subroutine ec_a



   !> Taken from http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
   subroutine INIT_RAND_SEED()
      use iso_fortran_env, only: int64
      integer, allocatable :: seed(:)
      integer :: i, n, un, istat, dt(8), pid
      integer(int64) :: t

      call random_seed(size = n)
      allocate(seed(n))
      ! First try if the OS provides a random number generator
      open(newunit=un, file="/dev/urandom", access="stream", &
      form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
         read(un) seed
         close(un)
      else
         ! Fallback to OR:ing the current time and pid. The PID is
         ! useful in case one launches multiple instances of the same
         ! program in parallel.
         call system_clock(t)
         if (t == 0) then
            call date_and_time(values=dt)
            t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
            + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
            + dt(3) * 24_int64 * 60 * 60 * 1000 &
            + dt(5) * 60 * 60 * 1000 &
            + dt(6) * 60 * 1000 + dt(7) * 1000 &
            + dt(8)
         end if
         pid = getpid()
         t = ieor(t, int(pid, kind(t)))
         do i = 1, n
            seed(i) = lcg(t)
         end do
      end if
      call random_seed(put=seed)
      contains
      ! This simple PRNG might not be good enough for real work, but is
      ! sufficient for seeding a better PRNG.
      function lcg(s)
         integer :: lcg
         integer(int64) :: s
         if (s == 0) then
            s = 104729
         else
            s = mod(s, 4294967296_int64)
         end if
         s = mod(s * 279470273_int64, 4294967291_int64)
         lcg = int(mod(s, int(huge(0), int64)), kind(0))
      end function lcg
   end subroutine INIT_RAND_SEED

end module data_block



!> Invert the small matrix
subroutine ec_INVERT_MATRIX_SHORT(lens,matrix,vector,solution_vector)
   use data_block, only: error
   implicit none
   !! I/O:
   INTEGER, intent(in)           :: lens
   DOUBLE PRECISION, intent(inout)  :: matrix(lens,lens)
   ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for
   ! condensed species, the pis and the delta log(n)
   DOUBLE PRECISION, intent(in)  :: vector(lens)
   double precision, intent(out) :: solution_vector(lens)
   !! Internal:
   INTEGER                      :: index(lens)
   DOUBLE PRECISION             :: indic

   interface
   subroutine ec_LUDCMP(a,indx,d)
      double precision, intent(inout)  :: a(:,:)
      integer, intent(out)             :: indx(:)
      double precision, intent(out)    :: d

      double precision                 :: vv(size(a,1))
      double precision, parameter      :: TINY = 1.0e-20
      integer                          :: j, n, imax
   end subroutine ec_LUDCMP
   subroutine ec_LUBKSB(a,indx,b)
      double precision, intent(in) :: a(:,:)
      integer, intent(in) :: indx(:)
      double precision, intent(inout) :: b(:)
      integer :: i,n,ii,ll
      double precision :: summ
   end subroutine ec_LUBKSB
   end interface

   solution_vector = vector

   call ec_LUDCMP(matrix,index,indic)
   if (error) RETURN
   call ec_LUBKSB(matrix,index,solution_vector)

end subroutine ec_INVERT_MATRIX_SHORT

!> Invert the big matrix
subroutine ec_INVERT_MATRIX_LONG(lens,matrix,vector,solution_vector)
   use data_block, only: N_gas, N_ions, remove_ions, reac_ion, error
   implicit none
   !! I/O:
   INTEGER, intent(in)           :: lens
   DOUBLE PRECISION, intent(inout)  :: matrix(lens,lens)
   double precision, intent(in)  :: vector(lens)
   double precision, intent(out) :: solution_vector(lens)
   ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for
   ! condensed species, the pis and the delta log(n)

   double precision              :: matrix_nions(lens-N_ions,lens-N_ions)
   DOUBLE PRECISION              :: vector_nions(lens-N_ions), &
   solution_vector_nions(lens-N_ions)

   !! Internal:
   INTEGER                       :: index(lens), corrf_i, corrf_j, index_nions(lens-N_ions)
   DOUBLE PRECISION              :: indic
   INTEGER                       :: i_mat, j_mat

   interface
   subroutine ec_LUDCMP(a,indx,d)
      double precision, intent(inout)  :: a(:,:)
      integer, intent(out)             :: indx(:)
      double precision, intent(out)    :: d

      double precision                 :: vv(size(a,1))
      double precision, parameter      :: TINY = 1.0e-20
      integer                          :: j, n, imax
   end subroutine ec_LUDCMP
   subroutine ec_LUBKSB(a,indx,b)
      double precision, intent(in) :: a(:,:)
      integer, intent(in) :: indx(:)
      double precision, intent(inout) :: b(:)
      integer :: i,n,ii,ll
      double precision :: summ
   end subroutine ec_LUBKSB
   end interface

   solution_vector = vector

   !matrix(N_reactants+N_atoms+1,N_reactants+N_atoms+1)
   IF (remove_ions) THEN
      vector_nions = 0d0
      matrix_nions = 0d0
      corrf_i = 0
      DO i_mat = 1, lens
         corrf_j = 0
         IF (i_mat <= N_gas) THEN
            IF (reac_ion(i_mat)) THEN
               corrf_i = corrf_i + 1
               cycle
            END IF
         END IF
         DO j_mat = 1, lens
            IF (j_mat <= N_gas) THEN
               IF (reac_ion(j_mat)) THEN
                  corrf_j = corrf_j + 1
                  cycle
               END IF
            END IF
            matrix_nions(j_mat-corrf_j,i_mat-corrf_i) = matrix(j_mat,i_mat)
         END DO
         vector_nions(i_mat-corrf_i) = vector(i_mat)
      END DO
      solution_vector_nions = vector_nions

      call ec_LUDCMP(matrix_nions,index_nions,indic)
      if (error) RETURN
      call ec_LUBKSB(matrix_nions,index_nions,solution_vector_nions)
      if (error) RETURN

      corrf_i = 0
      DO i_mat = 1, lens
         IF (i_mat <= N_gas) THEN
            IF (reac_ion(i_mat)) THEN
               corrf_i = corrf_i + 1
               cycle
            END IF
         END IF
         solution_vector(i_mat) = solution_vector_nions(i_mat-corrf_i)
      END DO
   ELSE
      call ec_LUDCMP(matrix,index,indic)
      if (error) RETURN
      call ec_LUBKSB(matrix,index,solution_vector)
      if (error) RETURN
   END IF

end subroutine ec_INVERT_MATRIX_LONG

!> LU decomposition, from numerical recipes
subroutine ec_LUDCMP(a,indx,d)
   use data_block, only: error, err_msg
   implicit none
   double precision, intent(inout)  :: a(:,:)
   integer, intent(out)             :: indx(:)
   double precision, intent(out)    :: d

   double precision                 :: vv(size(a,1))
   double precision, parameter      :: TINY = 1.0e-20
   integer                          :: j, n, imax, lu_asserteq

   n = lu_asserteq(size(a,1),size(a,2),size(indx),'ec_LUDCMP')
   if (error) RETURN

   d = 1.0
   vv = maxval(abs(a),dim=2)
   if (any(vv == 0.0)) then
      error = .true.
      err_msg = "ERROR: singular matrix in LUDCMP"
      RETURN
   end if
   vv = 1.0/vv
   do j = 1,n
      imax = (j-1)+maxloc(vv(j:n)*abs(a(j:n,j)), dim=1)
      if (j /= imax) then
         call lu_swap(n,a(imax,:),a(j,:))
         d = -d
         vv(imax) = vv(j)
      end if
      indx(j) = imax
      if (a(j,j) == 0.0) a(j,j) = TINY
      a(j+1:n,j) = a(j+1:n,j) / a(j,j)
      ! a(j+1:n,j+1:n) = a(j+1:n,j+1:n)-lu_outerprod(n-j,a(j+1:n,j),n-j,a(j,j+1:n))
      a(j+1:n,j+1:n) = a(j+1:n,j+1:n) - &
      (spread(a(j+1:n,j), dim=2, ncopies=n-j) * spread(a(j,j+1:n), dim=1, ncopies=n-j))
   end do
END SUBROUTINE ec_LUDCMP

!> LU back substitution
SUBROUTINE ec_LUBKSB(a,indx,b)
   use data_block, only: error
   implicit none
   double precision, intent(in) :: a(:,:)
   integer, intent(in) :: indx(:)
   double precision, intent(inout) :: b(:)
   integer :: i,n,ii,ll, lu_asserteq
   double precision :: summ

   n=lu_asserteq(size(a,1),size(a,2),size(indx),'ec_LUBKSB')
   if (error) RETURN

   ii=0
   do i=1,n
      ll=indx(i)
      summ=b(ll)
      b(ll)=b(i)
      if (ii /= 0) then
         summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
      else if (summ /= 0.0) then
         ii=i
      end if
      b(i)=summ
   end do
   do i=n,1,-1
      b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
   end do
END SUBROUTINE ec_LUBKSB

subroutine lu_swap(len,arr1,arr2)
   implicit none
   integer, intent(in)              :: len
   double precision, intent(inout)  :: arr1(len), arr2(len)
   double precision                 :: tamp(len)

   tamp = arr1
   arr1 = arr2
   arr2 = tamp
end subroutine lu_swap

function lu_asserteq(n1,n2,n3,label) result(m)
   use data_block, only: error, err_msg
   implicit none
   integer, intent(in)        :: n1, n2, n3
   character(len=9), intent(in)   :: label
   integer                    :: m

   if (n1==n2 .and. n2==n3) then
      m = n1
   else
      error = .true.
      err_msg = 'ERROR in '//trim(adjustl(label))//": the sizes of the matrix' don't add up..."
      ! print *, 'nrerror: assert_eq failed in: ', label
      ! STOP 'Program terminated by lu_asserteq'
      RETURN
   end if
end function lu_asserteq

! function lu_outerprod(len1,arr1,len2,arr2) result(arr3)
!    implicit none
!    integer, intent(in)           :: len1, len2
!    double precision, intent(in)  :: arr1(len1), arr2(len2)
!    double precision              :: arr3(len1, len2)

!    arr3 = spread(arr1, dim=2, ncopies=len2) * spread(arr2, dim=1, ncopies=len1)
! end function lu_outerprod