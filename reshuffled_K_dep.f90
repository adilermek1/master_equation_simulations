PROGRAM RK_Solution
	use barrier_energies_module
	use mpi
! Solve the differential equation dC/dt = ..  using rk4 subroutine
	implicit none
	
	integer, parameter :: barrier_rows=21
	integer, parameter :: barrier_columns = 65	
	real*8, dimension(barrier_rows, barrier_columns) :: barrier_energies_matrix
	character(len=200) :: barrier_directory, barrier_filename, barrier_filepath

	
	real*8 :: t0, t_final, t, h
	real*8, dimension(:), allocatable:: C, dCdt, Cout, kdecs_per_s, C_separate, dCdt_separate, C_equilibrium_separate, part_funcs, &
	eq_conc
	character(100) :: output, output1, output2, output3, output4, output5, output6, output7, output8, output9, output10, output11
	integer :: iunit, junit, kunit, nunit, munit, aunit, bunit, cunit, dunit, eunit, funit, j, unit, n, last_bound, k, gunit
	integer :: unit_K0, unit_K1, unit_K2, unit_K3, unit_K4, unit_K5, unit_K6, unit_K7, &
	unit_K8, unit_K9, unit_K10, unit_K11, unit_K12, unit_K13, unit_K14, unit_K15, unit_K16, &
	unit_K17, unit_K18, unit_K19, unit_K20, a, counter
	
! Parameters for kinetics
	character(len=100) :: directory, filename 
	
	character(len=200) :: filepath, filepath_K0, filepath_K1, filepath_K2, filepath_K3, filepath_K4, filepath_K5, filepath_K6, &
	filepath_K7, filepath_K8, filepath_K9, filepath_K10, filepath_K11, filepath_K12, filepath_K13, filepath_K14, filepath_K15, &
	filepath_K16, filepath_K17, filepath_K18, filepath_K19, filepath_K20
	
	integer :: num_states, i, ios, num_counter, num_limit, ios_K0, ios_K1, ios_K2, ios_K3, ios_K4, ios_K5, ios_K6, ios_K7, ios_K8, &
	ios_K9, ios_K10, ios_K11, ios_K12, ios_K13, ios_K14, ios_K15, ios_K16, ios_K17, ios_K18, ios_K19, ios_K20
	
	real*8, dimension(:), allocatable :: Energies, Gammas, Covalent_All, Vdw_All, Infinity_All, equilibrium_constants_m3, eq_concs
	real*8, dimension(:, :), allocatable :: truncated_matrix
	real*8, dimension(:, :), allocatable :: transition_matrix, Second_term_2d
	real*8, dimension(:, :, :), allocatable :: transition_matrix_3D
	integer, dimension(:, :), allocatable :: state_pointer
	real*8, dimension(:), allocatable :: First_term, Second_term, Second_term_truncated, sum_rows, &
	threshold_Energies_K, sum_rows_truncated, Second_term_main
	real*8 :: h_j_per_s, hbar_js, c_cm_per_s, C_O_per_m3, C_O2_per_m3, temp_k
	real*8 :: part_funcs_o2_per_m3
	character(len=3) :: o3_molecule
    character(len=2) :: o2_molecule, o_atom
	integer :: J_rot_start
	real*8 :: m_per_a0, j_per_k, j_per_cm, cm_per_k, kt_energy_cm, M_per_m3, ref_pressure_per_m3, pressure_ratio
	real*8 :: k0_m3_per_s, sigma0_m2, pi
	real*8, dimension(2) :: dE
	real*8 :: C_init, C_tot, dCdt_tot, K_eqs_tot, threshold_E
	real*8 :: Energy_tmp, Total_Gamma_tmp, Covalent_All_tmp, Vdw_All_tmp, Infinity_All_tmp
	real*8 :: k_rec, lower_energy_lim, upper_energy_lim, upper_gamma_lim, lower_gamma_lim
	integer :: threshold_j, Js, Ks, vib_sym_well, iteration_counter, print_freq, K_initial, K_final, band_width, K_exact
	integer, dimension(:), allocatable :: num_states_K, K_value, threshold_j_values, Resonance, num_counter_K, &
	unit_K, ios_K, istart, ifinish, cutoff_count, J_value, sym_value
	
	real*8 :: start_time, end_time, end_time1, end_time2, end_time3, tmp_energies, tmp_gammas, tmp_covalent_all, &
	tmp_vdw_all, tmp_infinity_all, tmp_k_value, tmp_resonance, kt_energy_j, dE_down, kij_min, start_time_matrix, &
	end_time_matrix, start_time_master_eq, end_time_master_eq, C_initial_O3_per_m3, &
	sum_eq_conc, E_dd, tmp_sym_value, tmp_j_value, w_1
	
	integer :: Ks_indep, cutoff_count_max, pe_num, ierr_a, myid, chunk_size, K_step, K_final_last, vib_sym_well_start, &
	vib_sym_well_last, J_initial, J_final, w
	
	logical :: print_detail, truncation, K_dependent, print_transition_matrix
	integer, dimension(:, :, :), allocatable :: num_states_J_K_sym, &
	num_counter_J_K_sym, unit_J_K_sym, ios_J_K_sym
	real*8, dimension(:, :, :), allocatable :: threshold_Energies_J_K_sym	
	
	character(len=200), dimension(:, :, :), allocatable :: filepath_JKsym
	character(len=10) :: Js_string, Ks_string
!---------------------------------------------------------------------------------------------------------------------------	
	call MPI_INIT(ierr_a)
	call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr_a)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, pe_num, ierr_a)
!	myid = 0
	if (myid == 0) print *, pe_num, "Number of processors"
		
	call cpu_time(start_time)
	
	barrier_directory = '/mmfs1/home/3436yermeka/kinetics/barrier_energies'
	barrier_filename = 'barrier_energies.csv'
	barrier_filepath = trim(barrier_directory) // '/' // trim(barrier_filename)
	
	barrier_energies_matrix = read_matrix(barrier_filepath)
	
	! Print the matrix
!	if (myid == 0) then
!	do i = 1, barrier_rows
!			write(*, '(65(1x,E21.14))') barrier_energies_matrix(i, :)
!	end do
!	end if
	
! Conversion factors	
	!m_per_a0 = get_m_per_a0();
	!j_per_k = get_j_per_k()
	!j_per_k = get_j_per_k()
	!j_per_cm = get_j_per_cm()
	m_per_a0 = 5.2917721092e-11
	j_per_k = 1.3806488e-23
	j_per_cm = 1.986445682933742e-23
	cm_per_k = j_per_k / j_per_cm
		
! Plank constant and speed of light c
	h_j_per_s = 6.62607015e-34
	hbar_js = 1.054571726e-34
	c_cm_per_s = 2.9979245800e+10
	pi = 4.D0*DATAN(1.D0)
	
! Gas Mixture Parameters for reagents and bath gas	
	C_O_per_m3 = 6.44e+18
	C_O2_per_m3 = 6.44e+20
	M_per_m3 = 10000*6.44e+24
	C_initial_O3_per_m3 = 0 !82510272735.4575
	
! Pressure related parameters
	ref_pressure_per_m3 = 6.44e+24
	pressure_ratio = M_per_m3 / ref_pressure_per_m3
	
	temp_k = 298
	sigma0_m2 = 2300 * m_per_a0**2

! Parameters for energy transfer model	
	dE_down = -43.13
	!dE(1) = -200.13 * j_per_cm
	!dE(2) = get_dE_up(dE(1), temp_k)
	
	j_per_k = get_j_per_k()
	kt_energy_j = temp_k * j_per_k;
	kt_energy_cm = temp_k * cm_per_k
	
! Energy Spectrum Parameters in wave numbers	
	lower_energy_lim = -30
	upper_energy_lim = 30
	upper_gamma_lim = 200
	lower_gamma_lim = 0.d0 ! 10**(-2)

! Ozone isotope Labels	
	o3_molecule = "666"
	o2_molecule = "66"
	o_atom = "6"

! Computing k0
	k0_m3_per_s = get_k0_2(o3_molecule, temp_k, sigma0_m2)

! Initial and final time and time step
	t0 = 0.0d0
! 0.01*1000e-9 / pressure_ratio = 1e-12 in case of 10000std of M_per_m3
	t_final = 1.4e-12 !135e-14 ! 1.35*0.01*1000e-9 / pressure_ratio
	h = 0.001*1e-13 !100e-20 ! 1*1e-9 / (pressure_ratio)
	band_width = 10
	print_freq = 1
	print_detail = .True.
	truncation = .False.
	K_dependent = .False.
	print_transition_matrix = .False.
	
! Specify the directory and file name  
	filename = "state_properties.fwc"

	num_states = 0	
	J_initial = 24
	J_final = 24
	K_exact = 0
	
	K_initial = 0
	K_step = 1
	
	K_final_last = 20
	vib_sym_well_start = 0
	vib_sym_well_last = 0
	
	allocate(num_states_J_K_sym(J_final+1, K_final_last+1, vib_sym_well_last+1))
	num_states_J_K_sym = 0
	allocate(threshold_Energies_J_K_sym(J_final+1, K_final_last+1, vib_sym_well_last+1))
	threshold_Energies_J_K_sym = 0
	allocate(filepath_JKsym(J_final+1, K_final_last+1, vib_sym_well_last+1))
	
	do vib_sym_well = vib_sym_well_start, vib_sym_well_last
		do Js = J_initial, J_final
			write(Js_string, '(I0)') Js
			
			if (Js < 20) then
				K_final = K_exact !Js
			else
				K_final = K_exact !K_final_last
			end if
			
			do Ks = K_initial, K_final, K_step
				write(Ks_string, '(I0)') Ks
				directory = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/"				
				if (vib_sym_well == 0) then
					! vib_sym_well = 0
					if (mod(Ks, 2) == 0) then
						! K is even
						call get_threshold_energy_K(o3_molecule, o2_molecule, Ks, threshold_E, threshold_j)
						threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1) = threshold_E
						directory = trim(directory) // "half_integers/J_" // trim(adjustl(Js_string)) // "/K_" &
						// trim(adjustl(Ks_string)) // "/symmetry_1"
						! Construct the full file path
						filepath = trim(directory) // '/' // trim(filename)
						filepath_JKsym(Js+1, Ks+1, vib_sym_well+1) = filepath
!						print*, filepath
							
							! Open and read the file with kinetics data
							open(newunit=unit, file=filepath, status='old', action='read', iostat=ios)
							! Skip the first line
							read(unit, *)
							! Count the number of data points
							do
								read(unit, *, iostat=ios) Energy_tmp, Total_Gamma_tmp, Covalent_All_tmp, Vdw_All_tmp, Infinity_All_tmp
								if (ios /= 0) exit
								if (Energy_tmp > max(threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1), barrier_energies_matrix(Ks+1, Js+1)) &
								+ lower_energy_lim .and. &
									Energy_tmp < max(threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1), barrier_energies_matrix(Ks+1, Js+1)) &
									+ upper_energy_lim &
									.and. Total_Gamma_tmp <= upper_gamma_lim) then
									num_states = num_states + 1
									num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1) = num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1) + 1
								end if
							end do
							if (myid == 0) then
							write(*, *) "vib_sym_well = 0, Number of states with K = ", Ks, "is: ", num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1)
							end if
							close(unit)	
						
					else
						! K is odd
						call get_threshold_energy_K(o3_molecule, o2_molecule, Ks, threshold_E, threshold_j)
						threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1) = threshold_E
						directory = trim(directory) // "J_" // trim(adjustl(Js_string)) // "/K_" &
						// trim(adjustl(Ks_string)) // "/symmetry_1"
						! Construct the full file path
						filepath = trim(directory) // '/' // trim(filename)
						filepath_JKsym(Js+1, Ks+1, vib_sym_well+1) = filepath
!						print*, filepath
						
							! Open and read the file with kinetics data
							open(newunit=unit, file=filepath, status='old', action='read', iostat=ios)
							! Skip the first line
							read(unit, *)
							! Count the number of data points
							do
								read(unit, *, iostat=ios) Energy_tmp, Total_Gamma_tmp, Covalent_All_tmp, Vdw_All_tmp, Infinity_All_tmp
								if (ios /= 0) exit
								if (Energy_tmp > max(threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1), barrier_energies_matrix(Ks+1, Js+1)) &
								+ lower_energy_lim .and. &
									Energy_tmp < max(threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1), barrier_energies_matrix(Ks+1, Js+1)) &
									+ upper_energy_lim &
									.and. Total_Gamma_tmp <= upper_gamma_lim) then
									num_states = num_states + 1
									num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1) = num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1) + 1									
								end if
							end do
							if (myid == 0) then
							write(*, *) "vib_sym_well = 0, Number of states with K = ", Ks, "is: ", num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1)
							end if
							close(unit)							
						
					end if
				else
					! vib_sym_well = 1
					if (mod(Ks, 2) == 0) then
						! K is even
						call get_threshold_energy_K(o3_molecule, o2_molecule, Ks, threshold_E, threshold_j)
						threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1) = threshold_E						
						directory = trim(directory) // "J_" // trim(adjustl(Js_string)) // "/K_" &
						// trim(adjustl(Ks_string)) // "/symmetry_1"
						! Construct the full file path
						filepath = trim(directory) // '/' // trim(filename)
						filepath_JKsym(Js+1, Ks+1, vib_sym_well+1) = filepath
!						print*, filepath

							! Open and read the file with kinetics data
							open(newunit=unit, file=filepath, status='old', action='read', iostat=ios)
							! Skip the first line
							read(unit, *)
							! Count the number of data points
							do
								read(unit, *, iostat=ios) Energy_tmp, Total_Gamma_tmp, Covalent_All_tmp, Vdw_All_tmp, Infinity_All_tmp
								if (ios /= 0) exit
								if (Energy_tmp > max(threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1), barrier_energies_matrix(Ks+1, Js+1)) &
								+ lower_energy_lim .and. &
									Energy_tmp < max(threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1), barrier_energies_matrix(Ks+1, Js+1)) &
									+ upper_energy_lim &
									.and. Total_Gamma_tmp <= upper_gamma_lim) then
									num_states = num_states + 1
									num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1) = num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1) + 1									
								end if
							end do
							if (myid == 0) then
							write(*, *) "vib_sym_well = 1, Number of states with K = ", Ks, "is: ", num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1)	
							end if
							close(unit)	

					else
						! K is odd
						call get_threshold_energy_K(o3_molecule, o2_molecule, Ks, threshold_E, threshold_j)
						threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1) = threshold_E						
						directory = trim(directory) // "half_integers/J_" // trim(adjustl(Js_string)) // "/K_" &
						// trim(adjustl(Ks_string)) // "/symmetry_1"
						! Construct the full file path
						filepath = trim(directory) // '/' // trim(filename)
						filepath_JKsym(Js+1, Ks+1, vib_sym_well+1) = filepath
!						print*, filepath
						
							! Open and read the file with kinetics data
							open(newunit=unit, file=filepath, status='old', action='read', iostat=ios)
							! Skip the first line
							read(unit, *)
							! Count the number of data points
							do
								read(unit, *, iostat=ios) Energy_tmp, Total_Gamma_tmp, Covalent_All_tmp, Vdw_All_tmp, Infinity_All_tmp
								if (ios /= 0) exit
								if (Energy_tmp > max(threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1), barrier_energies_matrix(Ks+1, Js+1)) &
								+ lower_energy_lim .and. &
									Energy_tmp < max(threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1), barrier_energies_matrix(Ks+1, Js+1)) &
									+ upper_energy_lim &
									.and. Total_Gamma_tmp <= upper_gamma_lim) then
									num_states = num_states + 1
									num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1) = num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1) + 1									
								end if
							end do
							if (myid == 0) then
							write(*, *) "vib_sym_well = 1, Number of states with K = ", Ks, "is: ", num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1)
							end if
							close(unit)
					end if
				end if
			end do
		end do
	end do
	
	if (myid == 0) then
	write(*, *) "Total number of states", num_states
	end if
	
! Allocate and store filtered data
	allocate(Energies(num_states))
	allocate(Gammas(num_states))
	allocate(Covalent_All(num_states))
	allocate(Vdw_All(num_states))
	allocate(Infinity_All(num_states))
	allocate(J_value(num_states))
	allocate(K_value(num_states))
	allocate(sym_value(num_states))
	allocate(Resonance(num_states))	
	
	Energies = 0
	Gammas = 0
	Covalent_All = 0
	Vdw_All = 0
	Infinity_All = 0
	J_value = 0
	K_value = 0
	sym_value = 0
	Resonance = 0
	
! Temporary variables for sorting	
	tmp_energies = 0
	tmp_gammas = 0
	tmp_covalent_all = 0
	tmp_vdw_all = 0
	tmp_infinity_all = 0
	tmp_k_value = 0
	tmp_resonance = 0	
	
	allocate(num_counter_J_K_sym(J_final+1, K_final_last+1, vib_sym_well_last+1))
	allocate(unit_J_K_sym(J_final+1, K_final_last+1, vib_sym_well_last+1))
	allocate(ios_J_K_sym(J_final+1, K_final_last+1, vib_sym_well_last+1))
	num_counter_J_K_sym = 1
	num_counter = 1
	
! Opening all files simultaneously	
	do vib_sym_well = vib_sym_well_start, vib_sym_well_last
	
		do Js = J_initial, J_final
		
			if (Js < 20) then
				K_final = K_exact !Js
			else
				K_final = K_exact !K_final_last
			end if
			
			do Ks = K_initial, K_final, K_step
				open(newunit=unit_J_K_sym(Js+1, Ks+1, vib_sym_well+1), file=filepath_JKsym(Js+1, Ks+1, vib_sym_well+1), &
				status='old', iostat=ios_J_K_sym(Js+1, Ks+1, vib_sym_well+1))
				read(unit_J_K_sym(Js+1, Ks+1, vib_sym_well+1), *)			
			end do
		end do
	end do
	
	
! Read and store filtered data
	do vib_sym_well = vib_sym_well_start, vib_sym_well_last
	
		do Js = J_initial, J_final
		
			if (Js < 20) then
				K_final = K_exact !Js
			else
				K_final = K_exact !K_final_last
			end if
			
! Read and store the filtered data
		do
		if (myid == 0) then
		print*, 'entering most outer do loop (sorting step)'
		end if
			do Ks = K_initial, K_final, K_step	
				do 
					if (num_counter_J_K_sym(Js+1, Ks+1, vib_sym_well+1) > num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1)) exit
					read(unit_J_K_sym(Js+1, Ks+1, vib_sym_well+1), *, iostat=ios_J_K_sym(Js+1, Ks+1, vib_sym_well+1)) &
					Energies(num_counter), Gammas(num_counter), &
					Covalent_All(num_counter), Vdw_All(num_counter), Infinity_All(num_counter)
					if (Energies(num_counter) > max(threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1), barrier_energies_matrix(Ks+1, Js+1)) &
					+ lower_energy_lim .and. &
					Energies(num_counter) < max(threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1), barrier_energies_matrix(Ks+1, Js+1)) &
						+ upper_energy_lim .and. Gammas(num_counter) <= upper_gamma_lim) then					
! Neglect Gamma if state is below threshold, otherwise call it a resonance					
						if (Energies(num_counter) < threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1)) then
							Gammas(num_counter) = 0
						else 
							Resonance(num_counter) = 1
						end if
!						print*, Energies(num_counter), Gammas(num_counter)							
!						Gammas(num_counter) = 0 ! Modification
						J_value(num_counter) = Js
						K_value(num_counter) =  Ks! Change this in the future
						sym_value(num_counter) = vib_sym_well
						if (num_counter > 1) then
							do k = num_counter, 2, -1
								if (Energies(k) < Energies(k-1)) then
									! Sorting Energies
									tmp_energies = Energies(k)
									Energies(k) = Energies(k-1)
									Energies(k-1) = tmp_energies
									! Sorting Gammas
									tmp_gammas = Gammas(k)
									Gammas(k) = Gammas(k-1)
									Gammas(k-1) = tmp_gammas
									! Sorting cov probs
									tmp_covalent_all = Covalent_All(k)
									Covalent_All(k) = Covalent_All(k-1)
									Covalent_All(k-1) = tmp_covalent_all
									! Soring vdw probs
									tmp_vdw_all = Vdw_All(k)
									Vdw_All(k) = Vdw_All(k-1)
									Vdw_All(k-1) = tmp_vdw_all
									! Sorting inf probs
									tmp_infinity_all = Infinity_All(k)
									Infinity_All(k) = Infinity_All(k-1)
									Infinity_All(k-1) = tmp_infinity_all
									! Sorting J_value array
									tmp_j_value = J_value(k)
									J_value(k) = J_value(k-1)
									J_value(k-1) = tmp_j_value
									! Sorting K_value array
									tmp_k_value = K_value(k)
									K_value(k) = K_value(k-1)
									K_value(k-1) = tmp_k_value
									! Sorting sym_value array
									tmp_sym_value = sym_value(k)
									sym_value(k) = sym_value(k-1)
									sym_value(k-1) = tmp_sym_value
									! Sorting Resonance array
									tmp_resonance = Resonance(k)
									Resonance(k) = Resonance(k-1)
									Resonance(k-1) = tmp_resonance	
								else
									exit
								end if						
							end do
						end if						
						num_counter_J_K_sym(Js+1, Ks+1, vib_sym_well+1) = num_counter_J_K_sym(Js+1, Ks+1, vib_sym_well+1) + 1
						num_counter = num_counter + 1
					end if
				end do
			end do
			if (num_counter > num_states) goto 10
			exit
 		end do
						
		end do
	end do
	
10	do vib_sym_well = vib_sym_well_start, vib_sym_well_last
		do Js = J_initial, J_final
			if (Js < 20) then
				K_final = K_exact !Js
			else
				K_final = K_exact !K_final_last
			end if	
			do Ks =  K_initial, K_final, K_step
				close(unit_J_K_sym(Js+1, Ks+1, vib_sym_well+1))
			end do
		end do
	end do
	
	if (myid == 0) write(*, *) "Counter over number of states:", num_counter-1
	call cpu_time(end_time2)
	if (myid == 0) write(*, *) "Time for the second reading and sorting:", end_time2-end_time1, "seconds"

! Chunk size in derivs
	if (pe_num == 1) then
		chunk_size = num_states
	else
		if (mod(num_states, pe_num) == 0) then
			chunk_size = (num_states / pe_num)
		else
			chunk_size = int(num_states / pe_num)+1
		end if
	
	end if
	
	if (myid == 0 ) then
		w = int(((pe_num)*chunk_size - num_states) / chunk_size)
		w_1 = real(((pe_num - w) * chunk_size - num_states), kind=8) / real(chunk_size, kind=8) * 100.0
		print*, 'Chunk size:', chunk_size
		print*, 'The number of not loaded processors:', w
		print*, 'Percent load of the last loaded processor:', w_1
	end if
	
	
! Check that sorting is correct	
	do i = 2, num_states
		if (Energies(i) < Energies(i-1)) then
			print *, 'error in sorting'
		end if
	end do

! Printing partition functions of O+O2 system in K blocks	
	allocate(part_funcs(K_final + 1))
	part_funcs = 0
	do Ks = K_initial, K_final, K_step
		if (K_dependent .eqv. .False.) then
			Ks_indep = 0
		else 
			Ks_indep = Ks
		end if
		call get_threshold_energy_K(o3_molecule, o2_molecule, Ks_indep, threshold_E, threshold_j)
		J_rot_start = threshold_j
		part_funcs(Ks + 1) = calc_part_func_O2_per_m3_total_mol(temp_k, o2_molecule, o_atom, J_rot_start, Ks)
	end do

! Loop over Ks to compute K-dependent Partition functions, Equilibrium constants and decay rates for all states
	allocate(equilibrium_constants_m3(num_states))
	allocate(kdecs_per_s(num_states))

	do  num_counter = 1, num_states
		Js = J_value(num_counter)
		Ks = K_value(num_counter)
		if (K_dependent .eqv. .False.) then
			Ks_indep = 0
		else
			Ks_indep = Ks
		end if
		call get_threshold_energy_K(o3_molecule, o2_molecule, Ks_indep, threshold_E, threshold_j)
		J_rot_start = threshold_j
		part_funcs_o2_per_m3 = calc_part_func_O2_per_m3_total_mol(temp_k, o2_molecule, o_atom, J_rot_start, Ks)

		equilibrium_constants_m3(num_counter) = calculate_formation_decay_equilibrium(Energies(num_counter), &
		temp_k, part_funcs_o2_per_m3, Js, Ks)
		kdecs_per_s(num_counter) = Gammas(num_counter) * j_per_cm / hbar_js
	end do
	
	if (myid == 0) then
	output = 'spectrum_info.txt'
	open(newunit=iunit, file=output, status='replace')
	do j = 1, num_states
!		if (K_value(j) == 0) then
		write(iunit, '((I8, 1x, E19.12, 1x, E19.12, 1x, I8, 1x, I8, 1x, I8, 1x, I8, 1x, E19.12, 1x, E19.12, 1x, E19.12, 1x, E19.12))') j, &
			Energies(j), Gammas(j), J_value(j), K_value(j), sym_value(j), Resonance(j), Covalent_All(j), &
			Vdw_All(j), Infinity_All(j), equilibrium_constants_m3(j)			
			flush(iunit)
!		end if
	end do
	close(iunit)
	end if

! Sum of equilibrium constants
	K_eqs_tot = 0
	do  i = 1, num_states
		K_eqs_tot = K_eqs_tot + equilibrium_constants_m3(i)
	end do

	allocate(eq_concs(num_states))
	do i = 1, num_states
		eq_concs(i) = equilibrium_constants_m3(i) * C_O_per_m3 * C_O2_per_m3
	end do	
	
	call cpu_time(start_time_matrix)
	! Compute and print transition matrix	
	allocate(transition_matrix(chunk_size, num_states))
	transition_matrix = calculate_transition_matrix_unitless(Energies, Covalent_All, Vdw_All, Infinity_All, dE_down)

! Scaling of the matrices
	transition_matrix = transition_matrix * k0_m3_per_s * M_per_m3

!	output11 = 'distribution_of_concentrations_over_kij_kji.txt'
!	open(newunit=gunit, file=output11, status='replace')
!	do i = 1, num_states
!		do j = 1, num_states
!			if (i /= j) then
!				write(gunit, '(I8, 1x, I8, 1x, E19.12, 1x, E19.12, 1x, E19.12, 1x, E19.12)') i, j, &
!				eq_concs(i)*transition_matrix(i, j), eq_concs(j)*transition_matrix(j, i), &
!				(eq_concs(i)*transition_matrix(i, j)-eq_concs(j)*transition_matrix(j, i))/&
!				((0.5)*(eq_concs(i)*transition_matrix(i, j)+eq_concs(j)*transition_matrix(j, i)))
!			end if
!		end do
!	end do
!	close(gunit)
	
	
!	if (myid == 0) then
!		if (print_transition_matrix .eqv. .True.) then
!			output2 = 'transition_matrix.csv'
!			open(newunit=kunit, file=output2, status='replace')
!			do j = 1, size(transition_matrix, 1)
!				write(kunit, '(1x,I8,a1)', advance='no') j, ','
!			end do
!				write (kunit, *)
!			
!			do i = 1, size(transition_matrix, 1)
!				do j = 1, size(transition_matrix, 2)
!				 if (j == 1) then
!					write(kunit, '(I8,a1,1x,E19.12,a1,1x)', advance='no') i, ',', transition_matrix(i, j), ','
!				 else 
!					write(kunit, '(E19.12,a1,1x)', advance='no') transition_matrix(i, j), ',' 
!				 end if
!				end do
!				write (kunit, *)
!			end do
!			close(kunit)
!		end if
!	end if

!! Counting the elements which are greater than cuttoff value in the transition matrix
!	output7 = 'cuttoff_count_numbers.out'
!	open(newunit=cunit, file=output7, status='replace')
!	output8 = 'truncated_matrix.csv'
!	open(newunit=dunit, file=output8, status='replace')
!	output9 = 'state_pointer.csv'
!	open(newunit=eunit, file=output9, status='replace')
!	
!	cutoff_count_max = 0
!	kij_min = 0.1
!	
!	allocate(istart(num_states))
!	allocate(ifinish(num_states))
!	allocate(cutoff_count(num_states))
!	allocate(state_pointer(num_states, num_states))
!	allocate(truncated_matrix(num_states, num_states))
!	ifinish = num_states
!	cutoff_count = 0
!	state_pointer = 0
!	truncated_matrix = 0
!	
!	do i = 1, num_states
!!	  istart(i) = 0
!		do j = 1, num_states
!			 if (transition_matrix(j, i) .gt. kij_min) then 
!			 cutoff_count(i) = cutoff_count(i) + 1
!!			 truncated_matrix(cutoff_count(i), i) = transition_matrix(j, i)
!			 truncated_matrix(j, i) = transition_matrix(j, i)
!			 truncated_matrix(i, j) = transition_matrix(i, j)
!			 state_pointer(cutoff_count(i), i) = j		 
!			 ifinish(i) = j	
!!			 if (istart(i) .eq. 0) istart(i) = j			 
!			 end if
!		end do
!				
!		if (cutoff_count(i) .gt. cutoff_count_max) cutoff_count_max = cutoff_count(i)
!	end do
!
!! Finding istart indexes	
!	do i = 1, num_states
!		istart(i) = 0
!		do j = 1, num_states
!			if (truncated_matrix(j,i) /= 0) then
!				if (istart(i) .eq. 0) istart(i) = j
!			end if
!		end do
!		istart(1) = 1
!		if (myid == 0) write(cunit, *) i, cutoff_count(i), istart(i), ifinish(i), ifinish(i) - istart(i) + 1, &
!		ifinish(i) - istart(i) + 1 - cutoff_count(i)
!	end do
!	if (myid == 0) write(cunit, *) cutoff_count_max
!	close(cunit)
!
!! Print truncation matrix
!	if (myid == 0) then
!		if (print_transition_matrix .eqv. .True.) then	
!			do j = 1, cutoff_count_max
!				write(dunit, '(1x,I8,a1)', advance='no') j, ','
!			end do
!				write (dunit, *)
!			
!			do i = 1, size(truncated_matrix, 1)
!				do j = 1, cutoff_count_max
!				 if (j == 1) then
!					write(dunit, '(I8,a1,1x,E19.12,a1,1x)', advance='no') i, ',', truncated_matrix(j, i), ','
!				 else 
!					write(dunit, '(E19.12,a1,1x)', advance='no') truncated_matrix(j, i), ',' 
!				 end if
!				end do
!					write (dunit, *)
!			end do
!			close(dunit)
!		end if
!	end if
!
!! Print state_pointer matrix
!	if (myid == 0) then
!		if (print_transition_matrix .eqv. .True.) then	
!			do j = 1, cutoff_count_max
!				write(eunit, '(1x,I8,a1)', advance='no') j, ','
!			end do
!				write (eunit, *)
!			
!			do i = 1, size(state_pointer, 1)
!				do j = 1, cutoff_count_max
!				 if (j == 1) then
!					write(eunit, '(I8,a1,1x,I8,a1,1x)', advance='no') i, ',', state_pointer(j, i), ','
!				 else 
!					write(eunit, '(I8,a1,1x)', advance='no') state_pointer(j, i), ',' 
!				 end if
!				end do
!					write (eunit, *)
!			end do
!			close(eunit)
!		end if
!	end if
	
	
	
	call cpu_time(end_time_matrix)
	if (myid == 0) write(*, *) "Time for calculating and writing matricies (transition, truncation and state pointer):", &
	end_time_matrix-start_time_matrix, "seconds"
	
	
	
! Allocate global transition matrix
	if (myid .eq. 0) then
		allocate(transition_matrix_3D(chunk_size, num_states, pe_num))
		transition_matrix_3D = 0
	end if
	
! Gather local transition matrix from all processes to global transition matrix
	call MPI_GATHER(transition_matrix, chunk_size*num_states, MPI_REAL8, &
                  transition_matrix_3D(:, :, myid+1), chunk_size*num_states, MPI_REAL8, 0, MPI_COMM_WORLD, ierr_a)
	call MPI_BARRIER( MPI_COMM_WORLD, ierr_a )				  

! Print global_matrix to a file on the master process
	if (myid == 0) then
	  call write_global_matrix_to_file(transition_matrix_3D, "output_global_matrix.csv", pe_num)
	end if			  
	
	if (myid .eq. 0) deallocate(transition_matrix_3D)
	
! Scaling of the matrices
!	transition_matrix = transition_matrix * k0_m3_per_s * M_per_m3
	truncated_matrix = truncated_matrix * k0_m3_per_s * M_per_m3

! Compute sum of rows in the transition matrix
	allocate(sum_rows(num_states))
	allocate(sum_rows_truncated(num_states))
	sum_rows = 0
	sum_rows_truncated = 0
!	do i = 1, num_states
!		do j = 1, num_states
!				sum_rows(i) = sum_rows(i) + transition_matrix(i, j)
!				sum_rows_truncated(i) = sum_rows_truncated(i) + truncated_matrix(i, j)
!		end do
!	end do

	do i = myid*chunk_size+1, myid*chunk_size+chunk_size
	   a = i - myid*chunk_size
	   do j = 1, num_states
			sum_rows(i) = sum_rows(i) + transition_matrix(a, j)
	   end do
	end do
	
!  Compute Constant terms for Master-equations
	allocate(First_term(num_states))
	allocate(Second_term(num_states))
	allocate(Second_term_truncated(num_states))
	do  i = 1, num_states
		First_term(i) = kdecs_per_s(i)*equilibrium_constants_m3(i)*C_O_per_m3*C_O2_per_m3
		Second_term(i) = kdecs_per_s(i) + sum_rows(i)
!		Second_term_truncated(i) = kdecs_per_s(i) + sum_rows_truncated(i)
	end do
	
	do i = myid*chunk_size+1, myid*chunk_size+chunk_size
		a = i - myid*chunk_size
		Second_term(i) = kdecs_per_s(i) + sum_rows(i)
!		write(5000+myid, *) i, a, Second_term(i)
	end do
		
! Setting up and printing the initial conditions for concentrations at t=0
	allocate(C(num_states))
	allocate(dCdt(num_states))
	
	t = t0
	C = C_initial_O3_per_m3
!	do i = 1, num_states
!		C(i) = equilibrium_constants_m3(i)*C_O_per_m3*C_O2_per_m3
!	end do

	if (myid == 0) then
	output1 = 'propagated_concentration.txt'
	open(newunit=junit, file=output1, status='replace')
	if (print_detail .eqv. .True.) then
		write(junit, '(*(E19.12,1x))') t, C
	end if
	end if
	
	iteration_counter = 0

! Compute Initial dCdt, Total dCdt and krec
	call derivs(t, C, dCdt)
!	if (myid == 0) write(*, '(*(E19.12,1x))', advance='no') dCdt

	C_tot = 0
	dCdt_tot = 0
		do  i = 1, num_states
		   C_tot = C_tot + C(i)
		   dCdt_tot = dCdt_tot + (Covalent_All(i) + Vdw_All(i))*dCdt(i)
!			if (myid == 0) write(3000, '(I8,1x,E20.12,1x,I8,1x,E20.12,1x,E20.12)')  myid+1, t, i, C(i), dCdt(i)
		end do
	k_rec = dCdt_tot / (M_per_m3*(C_O_per_m3*C_O2_per_m3 - C_tot/K_eqs_tot))
	
! Separate derivatives over K blocks at the initial moment of time 
	allocate(dCdt_separate(K_final+1))
	allocate(C_separate(K_final+1))
	allocate(C_equilibrium_separate(K_final+1))
	dCdt_separate = 0
	C_separate = 0
	C_equilibrium_separate = 0
	  	  
! Professor's method of separate sum over K blocks		 
	do i = 1, num_states
			dCdt_separate(K_value(i) + 1) = dCdt_separate(K_value(i) + 1) + (Covalent_All(i) + Vdw_All(i))*dCdt(i)
			C_separate(K_value(i) + 1) = C_separate(K_value(i) + 1) + C(i)
			C_equilibrium_separate(K_value(i) + 1) = C_equilibrium_separate(K_value(i) + 1) + &
				equilibrium_constants_m3(i) * C_O_per_m3 * C_O2_per_m3
	end do
		
	if (myid == 0) then
		output3 = 'recombination_coefficient_and_dCdt_tot.txt'
		open(newunit=nunit, file=output3, status='replace')
		write(nunit, '(*(E19.12,1x))') t, k_rec, dCdt_tot, dCdt_separate
			
		output4 = 'total_concentrations.txt'
		open(newunit=munit, file=output4, status='replace')
		write(munit, '(*(E19.12,1x))') t, C_tot, C_separate
		
		output5 = 'equilibrium_concentrations.txt'
		open(newunit=aunit, file=output5, status='replace')
		write(aunit, '(*(E19.12,1x))') h, 	 equilibrium_constants_m3 * C_O_per_m3 * C_O2_per_m3
		write(aunit, '(*(E19.12,1x))') t_final, equilibrium_constants_m3 * C_O_per_m3 * C_O2_per_m3
		close(aunit)
		
		output6 = 'equilibrium_concentrations_in_K_blocks'
		open(newunit=bunit, file=output6, status='replace')
		write(bunit, '(*(E19.12,1x))') h, 1.d0, C_equilibrium_separate
		write(bunit, '(*(E19.12,1x))') t_final, 1.d0, part_funcs
		close(bunit)
		
		output10 = 'equilibrium_concentrations_transitions_only.txt'
		open(newunit=funit, file=output10, status='replace')
		sum_eq_conc = 0
		do i = 1, num_states
			sum_eq_conc = sum_eq_conc + exp(-Energies(i)/kt_energy_cm)
		end do
		
		allocate(eq_conc(num_states))
		do i = 1, num_states
			eq_conc(i) = num_states*C_initial_O3_per_m3 * exp(-Energies(i)/kt_energy_cm) / sum_eq_conc
		end do
		write(funit, '(*(E19.12,1x))') h, 	 eq_conc
		write(funit, '(*(E19.12,1x))') t_final, eq_conc
		close(funit)		
		
		call cpu_time(end_time3)
		write(*, *) "All parameters before main propagation loop:", end_time3-start_time, "seconds"	
	end if
	
	call cpu_time(start_time_master_eq)
	
!  The main time-propogation loop
	allocate(Cout(num_states))
	do while (t<=t_final)
			  call derivs(t, C, dCdt)
			  
!			  if (iteration_counter == 1) then
!				  do i = 1, num_states
!!					 if (myid == 0) write(4000, '(I8,1x,E20.12,1x,I8,1x,E20.12,1x,E20.12,1x,E19.12)')  myid+1, t, i, C(i), dCdt(i), k_rec
!				  end do 
!			  end if		  
!			  if (iteration_counter == 20) stop

! Solve the differential equation using rk4 subroutine and update variables
			  call rk4(C, dCdt, num_states, t, h, Cout, derivs)
			    
			  t = t + h
			  C = Cout 	  
			  iteration_counter = iteration_counter + 1
			  			  
! Processing and printing of the results: Total dCdt and krec
			  if (mod(iteration_counter, print_freq) == 0) then			  
				  if (print_detail .eqv. .True.) then
					if (myid == 0) write(junit, '(*(E19.12,1x))') t, C
				  end if
				  
				  C_tot = 0
				  dCdt_tot = 0
				  do  i = 1, num_states
					C_tot = C_tot + C(i)
					dCdt_tot = dCdt_tot + (Covalent_All(i) + Vdw_All(i))*dCdt(i)
				  end do
				  
				  k_rec = dCdt_tot / (M_per_m3*(C_O_per_m3*C_O2_per_m3 - C_tot/K_eqs_tot))
! Separate derivatives over K blocks after propagation 	  
				  dCdt_separate = 0
			      C_separate = 0			 				  
! Professor's method of separate sum over K blocks
				  do i = 1, num_states
					 dCdt_separate(K_value(i) + 1) = dCdt_separate(K_value(i) + 1) + (Covalent_All(i) + Vdw_All(i))*dCdt(i)
					 C_separate(K_value(i) + 1) = C_separate(K_value(i) + 1) + C(i)
				  end do
				  
				  if (myid == 0) write(nunit, '(*(E19.12,1x))') t, k_rec, dCdt_tot, dCdt_separate
				  if (myid == 0) write(munit, '(*(E19.12,1x))') t, C_tot, C_separate
!				  flush(nunit)
!				  flush(munit)

			  end if
			  
	end do
	
	close(junit)
	close(nunit)
	close(munit)
	
	call cpu_time(end_time_master_eq)
	if (myid == 0) write(*, *) "Master equation propagation:", start_time_master_eq-end_time_master_eq, "seconds"
	
	call cpu_time(end_time)
	if (myid == 0) write(*, *) "CPU Time:", end_time-start_time, "seconds"

	call MPI_Finalize(ierr_a)
! Main part of the code ends here
!-------------------------------------------------------------------------------------------------------------

	contains
	
	subroutine write_global_matrix_to_file(matrix, filename, num_procs)
	  real*8, dimension(:, :, :) :: matrix
	  character(len=*), intent(in) :: filename
	  integer, intent(in) :: num_procs
	  integer :: unit, i, j, k, chunk_size_last

	  ! Open the file for writing
	  open(newunit=unit, file=filename, status='replace')

	  ! Print header
	  do j = 1, size(matrix, 2)
		write(unit, '(1x,I8,a1)', advance='no') j, ','
	  end do
	  write(unit, *)

	  ! Loop over each processor
	  do k = 1, num_procs
		! Loop over each row of the matrix and write to the file
		do i = (k-1)*chunk_size+1, (k-1)*chunk_size+chunk_size
		a = i - (k-1)*chunk_size
		
		if (myid .eq. pe_num-1) then
			chunk_size_last = num_states - (pe_num-1)*chunk_size
			if (a > chunk_size_last) exit
		end if
		
		if (i .gt. num_states) exit
		
		  ! Loop over each column of the matrix and write to the file
		  do j = 1, size(matrix, 2)
			if (j == 1) then
			  write(unit, '(I8,a1,1x,E19.12,a1,1x)', advance='no') i, ',', matrix(a, j, k), ','
			else
			  write(unit, '(E19.12,a1,1x)', advance='no') matrix(a, j, k), ','
			end if
		  end do
		  write(unit, *)
		end do
	  end do

	  ! Close the file
	  close(unit)
	end subroutine write_global_matrix_to_file
		 	
	subroutine derivs(t, C, dCdt)
		use mpi
		integer :: myid_counter, counter
		real*8, intent(in) :: t, C(num_states)
		real*8, intent(out) :: dCdt(num_states)
		real*8 :: dCdt_myid(chunk_size), dCdt_myid2(chunk_size, pe_num)
		real*8 Third_term
		integer :: chunk_size_last
		
		myid_counter = 0
		dCdt_myid = 0
! Calculate the derivatives dC/dt
! Loop over final states, same as equation numbers:
!		do i = myid+1, num_states, pe_num
		do i = myid*chunk_size+1, myid*chunk_size+chunk_size
		a = i - myid*chunk_size
		
		if (myid .eq. pe_num-1) then
			chunk_size_last = num_states - (pe_num-1)*chunk_size
			if (a > chunk_size_last) exit
!			if (a > chunk_size_last) goto 30
		end if
		
		myid_counter = myid_counter + 1
!		print*, a, myid_counter
!		do i = 1, num_states
! calculate sum for individual state - Third_term		
			Third_term = 0 ! Initialize Third_term for this state			
! Loop over initial states, all terms of the summ:
				if (truncation .eqv. .False.) then
					do j = 1, num_states 					
!						Third_term = Third_term + transition_matrix(j, a) * C(j)
						Third_term = Third_term + transition_matrix(a, j)*(equilibrium_constants_m3(i)/equilibrium_constants_m3(j))*C(j)
!					if (myid .eq. pe_num-1) stop 
!					write(myid*1000+1, '(I8,1x,I8,1x,E20.12,1x,E20.12,1x,I8)') i, j, transition_matrix(a, j)*(equilibrium_constants_m3(i)/equilibrium_constants_m3(j)), C(j), a	
					end do
									
				else
!					istart(i) = max(1, i-band_width/2)
!					ifinish(i) = min(i+band_width/2, num_states)
					do j = istart(i), ifinish(i) 
!					do j = 1, cutoff_count(i)
!						Third_term = Third_term + truncated_matrix(j, i) * C(state_pointer(j,i))
						Third_term = Third_term + transition_matrix(j, i) * C(j)
					end do
				end if
				
			if (truncation .eqv. .False.) then
				dCdt_myid(myid_counter) = First_term(i) - Second_term(i)*C(i) + Third_term
!				dCdt(i) = First_term(i) - Second_term(i)*C(i) + Third_term
!			write(myid*1000+1, '(I8,1x,E20.12,1x,E20.12,1x,E20.12,1x,E20.12,1x,I8,1x,E20.12)') i, First_term(i), Second_term(i), C(i), dCdt_myid(myid_counter), myid_counter, Third_term				
			else
				dCdt_myid(myid_counter) = First_term(i) - Second_term_truncated(i)*C(i) + Third_term			
!				dCdt(i) = First_term(i) - Second_term_truncated(i)*C(i) + Third_term
			end if
			
		end do
		
			call MPI_GATHER(dCdt_myid, chunk_size, MPI_REAL8, dCdt_myid2, chunk_size, MPI_REAL8, 0, MPI_COMM_WORLD, ierr_a)			
			call MPI_BARRIER( MPI_COMM_WORLD, ierr_a )
						
!			if (myid == 0) write(1000, *)  dCdt_myid2
			 
			if (myid == 0 ) then 
				do j = 1, pe_num
					do i = 1, chunk_size
						counter = (j-1)*chunk_size + i
						dCdt(counter) = dCdt_myid2(i, j)
!						write(1000, *)  dCdt(counter), counter						
					end do
				 end do
			end if 
					
			call MPI_BCAST(dCdt, num_states, MPI_REAL8, 0, MPI_COMM_WORLD, ierr_a)
			call MPI_BARRIER(MPI_COMM_WORLD, ierr_a)
						
!			if (myid == 0) write(1002, '(*(E19.12,1x))', advance='no') dCdt
!			dCdt = First_term - kdecs_per_s*C + MATMUL(modified_transition_matrix_per_s_transposed, C)
			
	end subroutine derivs
	
	subroutine rk4(y,dydx,n,x,h,yout,derivs)
! Taken from Numerical Recipes
	   integer n,NMAX  
	   real*8 h,x,dydx(n),y(n),yout(n)  
	   EXTERNAL derivs  
	   PARAMETER (NMAX=1000000)
	   integer i  
	   real*8 h6, hh, xh, dym(NMAX), dyt(NMAX), yt(NMAX)

	   hh=h*0.5d0  
	   h6=h/6.d0  
	   xh=x+hh  

	   do 11 i=1,n  
		 yt(i)=y(i)+hh*dydx(i)  
	11    continue  

	   call derivs(xh,yt,dyt)  
	   
	   do 12 i=1,n  
		 yt(i)=y(i)+hh*dyt(i)  
	12    continue  

	   call derivs(xh,yt,dym)  
	   
	   do 13 i=1,n  
		 yt(i)=y(i)+h*dym(i)  
		 dym(i)=dyt(i)+dym(i)  
	13    continue  

	   call derivs(x+h,yt,dyt)  

	   do 14 i=1,n  
		 yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.d0*dym(i))  
	14    continue  
	   return  
	end subroutine rk4
	  
	function calculate_transition_matrix_unitless(Energies, Cov, Vdw, Inf, dE_down) result(matrix)
! Calculates unitless state-to-state transition matrix (matrix(i, j) = kappa i->j)
! dE = [down, up]
		real*8, dimension(:), intent(in) :: Energies(:)
		real*8, dimension(:), intent(in) :: Cov(:), Vdw(:), Inf(:)
		real*8 :: dE_down, dE_up
		real*8 :: ptran, matrix_tmp
		real*8, allocatable, dimension(:, :) :: matrix
		integer :: i, j, a, chunk_size_last
	
!		allocate(matrix(num_states, num_states))
		allocate(matrix(chunk_size, num_states))
		matrix = 0
		
!		do j = 1, size(matrix, 2)
		do j = myid*chunk_size+1, myid*chunk_size+chunk_size
		a = j - myid*chunk_size
		
		if (j > num_states) exit
		
		if (myid .eq. pe_num-1) then
			chunk_size_last = num_states - (pe_num-1)*chunk_size
			if (a > chunk_size_last) goto 40
		end if		
		
		do i = 1, num_states
!		   do i = 1, j-1
			  ptran = Cov(i)*Cov(j) + Vdw(i)*Vdw(j)! + Inf(i)*Inf(j)
			  dE_up = get_dE_up(dE_down, temp_k, Energies(i), Energies(j), equilibrium_constants_m3(i), equilibrium_constants_m3(j))
			  if (Energies(j) > Energies(i)) then
				matrix(a, i) = ptran * exp((Energies(j) - Energies(i)) / dE_down)
			  else if (Energies(j) < Energies(i)) then
				matrix_tmp = ptran * exp((Energies(i) - Energies(j)) / dE_down)
				matrix(a, i) = matrix_tmp * eq_concs(i) / eq_concs(j)
!				matrix(a, i) = ptran * exp((Energies(j) - Energies(i)) / dE_up)
			  else 
				matrix(a, i) = 0
			  end if
		   end do
		end do

!40		call MPI_BARRIER( MPI_COMM_WORLD, ierr_a )

!! 		  Print header
!		  do i = 1, size(matrix, 2)
!			write(myid+100, '(1x,I8,a1)', advance='no') i, ','
!		  end do
!		  write(myid+100, *)
!
!			! Loop over each row of the matrix and write to the file
!			do j = myid*chunk_size+1, myid*chunk_size+chunk_size
!			a = j - myid*chunk_size
!			
!			if (myid .eq. pe_num-1) then
!				chunk_size_last = num_states - (pe_num-1)*chunk_size
!				if (a > chunk_size_last) exit
!			end if
!		
!			  ! Loop over each column of the matrix and write to the file
!			  do i = 1, size(matrix, 2)				
!				if (i == 1) then
!				  write(myid+100, '(I8,a1,1x,E19.12,a1,1x)', advance='no') j, ',', matrix(a, i), ','
!				else
!				  write(myid+100, '(E19.12,a1,1x)', advance='no') matrix(a, i), ','
!				end if
!			  end do
!			  write(myid+100, *)
!			end do
		
40	end function calculate_transition_matrix_unitless
	   
	subroutine readFWCFile(filepath, Energy_per_cm, Total_Gamma_per_cm, Covalent_All, Vdw_all, Infty, A_per_cm, &
			   B_per_cm, C_per_cm, numDataPoints)
		character(len=*), intent(in) :: filepath
		real*8, dimension(:), allocatable, intent(out) :: Energy_per_cm, Total_Gamma_per_cm
		real*8, dimension(:), allocatable, intent(out) :: Covalent_All, Vdw_all, Infty, A_per_cm, B_per_cm, C_per_cm  
		integer, intent(out) :: numDataPoints
		integer :: unit, i
		character(len=100) :: header
! open the file for reading
		open(newunit=unit, file=filepath, status='old', action='read')
! Skip the header line
		read(unit, *)
! Determine the number of data points
		numDataPoints = 0
		
		do 
		  read(unit, *, iostat=i)
		  IF(i /= 0) EXIT
		  numDataPoints = numDataPoints + 1
		end do
		
		! allocate memory for the PARAMETER values
		allocate(Energy_per_cm(numDataPoints))
		allocate(Total_Gamma_per_cm(numDataPoints))
		allocate(Covalent_All(numDataPoints))
		allocate(Vdw_all(numDataPoints))
		allocate(Infty(numDataPoints))
		allocate(A_per_cm(numDataPoints))
		allocate(B_per_cm(numDataPoints))
		allocate(C_per_cm(numDataPoints))
		! rewind the file
		rewind(unit)
		! Skip the header line again
		read(unit, *)
		! read the data into the array
		do i = 1, numDataPoints
		  read(unit, *) Energy_per_cm(i), Total_Gamma_per_cm(i), &
					  Covalent_All(i), Vdw_all(i), Infty(i), & 
					  A_per_cm(i), B_per_cm(i), C_per_cm(i) 
		end do
! close the file
		close(unit)
	end subroutine readFWCFile
	
!!!!!!!!!!!!!!!!!!!!!!! This block need to be finished
	function get_higher_barrier_threshold() result(threshold_energy_j)
		real*8 :: threshold_energy_j
		threshold_energy_j = 5.743224250413328e-23
	end function get_higher_barrier_threshold
	
	function get_threshold_energies_2() result(threshold_energy_j)
		real*8 :: threshold_energy_j
		threshold_energy_j = 5.743224250413328e-23
	end function get_threshold_energies_2
!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine get_threshold_energy_K(o3_molecule, o2_molecule, K, threshold_E, threshold_j) 
		character(len=3), intent(in) :: o3_molecule
		character(len=2), intent(in) ::  o2_molecule
		integer, intent(in) :: K
		real*8 :: threshold_energy_j
		real*8, intent(out) :: threshold_E
		integer, intent(out) :: threshold_j
		!real*8 :: channel_shift_j
		!channel_shift_j = get_channel_shift(o3_molecule, o2_molecule)
		threshold_j = get_o2_threshold_j(o2_molecule, K)
		threshold_energy_j = rigid_rotor_energy(threshold_j, get_inertia_moment_2(o2_molecule)) !+ channel_shift_j
		threshold_E = threshold_energy_j / j_per_cm
	end subroutine get_threshold_energy_K
	
	function get_o2_threshold_j(o2_molecule, K) result(o2_threshold_j)
		character(len=2)  o2_molecule
		integer :: K
		integer :: o2_threshold_j
		
		if (mod(K, 2) == 0 .and. is_monoisotopic(o2_molecule)) then
		  o2_threshold_j = K + 1
		else 
		  o2_threshold_j = K
		end if
		
	end function get_o2_threshold_j
	
	function get_inertia_moment_2(o2_molecule) result(I_kg_m2)
		character(len=2), intent(in) :: o2_molecule
		real*8 :: mu_rot_kg, I_kg_m2
		mu_rot_kg = get_mu_rot(o2_molecule)
		I_kg_m2 = get_inertia_moment(mu_rot_kg)
	end function get_inertia_moment_2
	
	function calculate_formation_decay_equilibrium(Energy, temp_k, part_funcs_o2_per_m3, J_value, K_value) result(Kfds_m3)
		real*8 :: temp_k, part_funcs_o2_per_m3, threshold_E, threshold_energy_j
		real*8 :: Energy, Total_Energy_j, J, j_per_k, kt_energy_j, Kfds_m3
		integer :: J_value, K_value
		
		call get_threshold_energy_K(o3_molecule, o2_molecule, K_value, threshold_E, threshold_j)
		j_per_k = get_j_per_k()
		kt_energy_j = temp_k * j_per_k;
		Total_Energy_j = Energy * j_per_cm
		threshold_energy_j = threshold_E * j_per_cm
			
		Kfds_m3 = (2*J_value+1)*exp(-(Total_Energy_j) / kt_energy_j) / part_funcs_o2_per_m3
		
	end function calculate_formation_decay_equilibrium
	
	function get_dE_up(dE_down, temp_k, Energy_i, Energy_j, K_eq_i, K_eq_j) result(dE_up)
! Calculates dE_up corresponding to dE_down to satisfy the reversibility principle
		real*8, intent(in) :: dE_down, temp_k, Energy_i, Energy_j, K_eq_i, K_eq_j
		real*8 :: j_per_k, kt_energy_j, dE_up , dE_down_j

		j_per_k = get_j_per_k()
		kt_energy_j = temp_k * j_per_k
		
!		dE_down_j = dE_down * j_per_cm
!		dE_up = (dE_down_j/j_per_cm) / (dE_down_j / kt_energy_j - 1)
		dE_up = 1.d0 / ( log(K_eq_j/K_eq_i)/(Energy_i-Energy_j) - 1.d0/dE_down )
		
	end function get_dE_up
	
	function calc_part_func_O2_per_m3_elec(temp_k) result(pfunc)
! returns electronic partition function of O2 (+O?) system
! The expression is taken from:
! Hathorn, B. C.; Marcus, R. A. 
! An Intramolecular Theory of the Mass-Independent Isotope Effect for Ozone. 
! II. Numerical Implementation at Low Pressures Using a Loose Transition State. 
! J. Chem. Phys. 2000, 113 (21), 9497–9509. 
! https://doi.org/10.1063/1.480267
		real*8, intent(in) :: temp_k
		real*8 :: pfunc
		real*8 :: j_per_k, j_per_cm, cm_per_k, kt_energy_cm
		j_per_k = get_j_per_k()
		j_per_cm = get_j_per_cm()
		cm_per_k = j_per_k / j_per_cm
		kt_energy_cm = temp_k * cm_per_k
		pfunc = 15 + 9 * exp(-158.5 / kt_energy_cm) + 3 * exp(-226.5 / kt_energy_cm)
	end function calc_part_funC_O2_per_m3_elec
	
	function get_j_per_k() result(j_per_k)
		real*8 :: j_per_k
		j_per_k = 1.3806488e-23
		end function get_j_per_k
		
		function get_j_per_cm() result(j_per_cm)
		real*8 :: j_per_cm
		j_per_cm = 1.986445682933742e-23
		end function get_j_per_cm
		
		function calc_part_func_O2_per_m3_vib(zpe_j, kt_energy_j) result(pfunc)
		real*8 :: zpe_j, kt_energy_j, pfunc
		pfunc = 1.0 / (1.0 - exp(-2.0 * zpe_j / kt_energy_j))
	end function calc_part_func_O2_per_m3_vib
	
	function calc_part_func_O2_per_m3_rot(mu_rot_kg, J_start, J_step, kt_energy_j, K) result(pfunc)
! mu_rot_kg - reduced mass of the O2 system
		integer, intent(in) :: J_start, J_step
		real*8, intent(in) :: mu_rot_kg, kt_energy_j
		real*8 :: eps, I_kg_m2, threshold_energy_j, energy_j, new_pfunc, pfunc, threshold_E
		integer :: J, threshold_j, K
		eps = 1e-10
		I_kg_m2 = get_inertia_moment(mu_rot_kg)
!		threshold_energy_j = rigid_rotor_energy(J_start, I_kg_m2)
		call get_threshold_energy_K(o3_molecule, o2_molecule, K, threshold_E, threshold_j)
		threshold_energy_j = threshold_E * j_per_cm
		
		pfunc = 0.0
		J = J_start
		do while (.true.)
		energy_j = rigid_rotor_energy(J, I_kg_m2)
		new_pfunc = pfunc + (2*J + 1) * exp(-(energy_j) / kt_energy_j)
		if (new_pfunc - pfunc < eps) then
		  exit
		end if
		pfunc = new_pfunc
		J = J + J_step
		end do
	end function calc_part_func_O2_per_m3_rot
	
	function get_inertia_moment(mu_rot_kg) result(I_kg_m2)
! Returns moment of inertia of a given O2 molecule
		real*8, intent(in) :: mu_rot_kg
		real*8 :: r0_m, I_kg_m2
		r0_m = get_o2_distance_a0() * get_m_per_a0()
		I_kg_m2 = mu_rot_kg * r0_m**2
	end function get_inertia_moment
	
	function get_o2_distance_a0() result(dist_a0)
! Returns equilibrium distance between atoms in an O2 molecule
		real*8 :: dist_a0
		dist_a0 = 2.2819
	end function get_o2_distance_a0
	
	function get_m_per_a0() result(m_per_a0)
! Returns number of meters in one Bohr
		real*8 :: m_per_a0
		m_per_a0 = 5.2917721092e-11
	end function get_m_per_a0

	function rigid_rotor_energy(J, I_kg_m2) result(energy_J)
! returns energy (in J) of a rigid rotor corresponding to given J and I (moment of inertia)
		integer J
		real*8 :: I_kg_m2, energy_J
		real*8 :: hbar_js
		hbar_js = get_hbar_js()
		energy_J = J * (J + 1) / (2 * I_kg_m2) * hbar_js**2
		end function rigid_rotor_energy
		
		function get_hbar_js() result(hbar_js)
		real*8 :: hbar_js
		hbar_js = 1.054571726d-34
	end function get_hbar_js
	
	function calc_part_func_O2_per_m3_trans(mu_trans_kg, kt_energy_j) result(pfunc)
! mu_trans_kg - reduced mass of the O+O2 system
		real*8 :: mu_trans_kg, kt_energy_j, hbar_js, de_broglie_wl, pfunc
		hbar_js = get_hbar_js()
		de_broglie_wl = sqrt(2 * pi / (mu_trans_kg * kt_energy_j)) * hbar_js
		pfunc = de_broglie_wl**(-3)
	end function calc_part_func_O2_per_m3_trans
						
	function calc_part_func_O2_per_m3_total(temp_k, zpe_j, mu_rot_kg, J_rot_start, J_rot_step, mu_trans_kg, K) result(pfunc)
! Calculates the overall partition function of O+O2 system
		real*8 :: temp_k, zpe_j, mu_rot_kg, mu_trans_kg
		real*8 :: j_per_k, kt_energy_j, pfunc_elec, pfunc_vib, pfunc_rot, pfunc_trans, pfunc
		integer J_rot_start, J_rot_step, K
		j_per_k = get_j_per_k()
		kt_energy_j = temp_k * j_per_k
		pfunc_elec = calc_part_func_O2_per_m3_elec(temp_k)
		pfunc_vib = calc_part_func_O2_per_m3_vib(zpe_j, kt_energy_j)
		pfunc_rot = calc_part_func_O2_per_m3_rot(mu_rot_kg, J_rot_start, J_rot_step, kt_energy_j, K)
		pfunc_trans = calc_part_func_O2_per_m3_trans(mu_trans_kg, kt_energy_j)
		pfunc = pfunc_elec * pfunc_vib * pfunc_rot * pfunc_trans
	end function calc_part_func_O2_per_m3_total
	
	function calc_part_func_O2_per_m3_total_mol(temp_k, o2_molecule, o_atom, J_rot_start, K) result(pfunc_m_3)
		real*8 :: temp_k, pfunc_m_3
		character(len=*), intent(in) :: o2_molecule, o_atom
		real*8 :: zpe_j, mu_rot_kg, mu_trans_kg
		integer J_rot_start, J_rot_step, K
		
		zpe_j = get_channel_zpe(o2_molecule)
		mu_rot_kg = get_mu_rot(o2_molecule)
		mu_trans_kg = get_mu_trans(o2_molecule, o_atom)
		if (is_monoisotopic(o2_molecule)) then
			J_rot_step = 2
			else
			J_rot_step = 1
		end if
		
		pfunc_m_3 = calc_part_func_O2_per_m3_total(temp_k, zpe_j, mu_rot_kg, J_rot_start, J_rot_step, mu_trans_kg, K)
	end function calc_part_func_O2_per_m3_total_mol
	
	function is_monoisotopic(molecule) result(res)
	  character(len=*), intent(in) :: molecule
	  logical :: res
	  integer :: i, molecule_len

	  molecule_len = len(molecule)
	  res = .true.
	  do i = 2, molecule_len
		if (molecule(i:i) /= molecule(1:1)) then
		  res = .false.
		  exit
		end if
	  end do
	end function is_monoisotopic
	
	function get_channel_zpe(o2_molecule) result(zpe_j)
! Returns zero-point energy of a given channel
		character(len=*), intent(in) :: o2_molecule
		real*8 :: zpe_cm, zpe_j

		if (o2_molecule == "66") then
		zpe_cm = 7.916382691754641e+02
		elseif (o2_molecule == "67") then
		zpe_cm = 7.799050607081324e+02
		elseif (o2_molecule == "68") then
		zpe_cm = 7.693708543295361e+02
		elseif (o2_molecule == "77") then
		zpe_cm = 7.679904668774815e+02
		elseif (o2_molecule == "78") then
		zpe_cm = 7.572885872707559e+02
		elseif (o2_molecule == "88") then
		zpe_cm = 7.464315071358510e+02
		endif

		zpe_j = zpe_cm * get_j_per_cm()
	end function get_channel_zpe
			
	function get_mu_rot(o2_molecule) result(mu_rot_kg)
! Returns reduced mass of a given O2 molecule
		character(len=*), intent(in) :: o2_molecule
		real*8 :: atom_masses_kg(2), mu_rot_kg

		atom_masses_kg = get_atom_masses(o2_molecule)
		mu_rot_kg = product(atom_masses_kg) / sum(atom_masses_kg)
	end function get_mu_rot
			
	function get_atom_masses(molecule) result(atom_masses)
		character(len=*), intent(in) :: molecule
		real*8, dimension(:), allocatable :: atom_masses
		real*8, dimension(3) :: oxygen_kg
		integer :: i, len

		oxygen_kg = get_oxygen_mass_amu() * get_kg_per_amu()
		len = len_trim(molecule)
		allocate(atom_masses(len))

		do i = 1, len
		select case (molecule(i:i))
		case ('6')
		  atom_masses(i) = oxygen_kg(1)
		case ('7')
		  atom_masses(i) = oxygen_kg(2)
		case ('8')
		  atom_masses(i) = oxygen_kg(3)
		end select
		end do
	end function get_atom_masses
				
	function get_oxygen_mass_amu() result(oxygen_amu)
! Returns masses of oxygen isotopoes (in amu)
		real*8, dimension(3) :: oxygen_amu

		oxygen_amu = [15.99491461956, 16.9991317, 17.9991596129]
	end function get_oxygen_mass_amu
	
	function get_nitrogen_mass_amu() result(nitrogen_amu)
! Returns masses of nitrogen isotopoes (in amu)
		real*8,dimension(2) :: nitrogen_amu
		
		nitrogen_amu = [14.0030740048, 15.0001088982]
		end function get_nitrogen_mass_amu
						
		function get_kg_per_amu() result(kg_per_amu)
		real*8 :: kg_per_amu

		kg_per_amu = 1.660538921d-27
		end function get_kg_per_amu
		
		function get_k0_2(o3_molecule, temp_k, sigma0_m2) result(k0_m3_per_s)
		real*8, intent(in) :: temp_k, sigma0_m2
		character(len=3), intent(in) :: o3_molecule
		real*8 :: ozone_mass_kg, third_body_mass_kg, kt_energy_j, k0_m3_per_s, kg_per_amu
		real*8, dimension(2) :: nitrogen_amu
		
		j_per_k = get_j_per_k();
		kg_per_amu = get_kg_per_amu();
		nitrogen_amu = get_nitrogen_mass_amu();
		
		ozone_mass_kg = 7.96805347915539e-26
		third_body_mass_kg = nitrogen_amu(1) * 2 * kg_per_amu
		kt_energy_j = temp_k * j_per_k
		k0_m3_per_s = get_k0(ozone_mass_kg, third_body_mass_kg, kt_energy_j, sigma0_m2)
		
	end function get_k0_2
	
	function get_k0(ozone_mass_kg, third_body_mass_kg, kt_energy_j, sigma_stab_m2) result(k_stab_m3_per_s)
! Calculates kstab corresponding to sigma_stab (stabilization cross-section)
! ozone_mass is the sum of masses of individual atoms, NOT the reduced mass (mu)
		real(8), intent(in) :: ozone_mass_kg, third_body_mass_kg, kt_energy_j, sigma_stab_m2
		real(8) :: mu_stab_kg, velocity_m_per_s, k_stab_m3_per_s

		mu_stab_kg = ozone_mass_kg * third_body_mass_kg / (ozone_mass_kg + third_body_mass_kg)
		velocity_m_per_s = sqrt(8 * kt_energy_j / (pi * mu_stab_kg))
		k_stab_m3_per_s = sigma_stab_m2 * velocity_m_per_s
		end function get_k0
							
		function get_mu_trans(o2_molecule, o_atom) result(mu_trans_kg)
		character(len=*), intent(in) :: o2_molecule, o_atom
		real*8 :: mu_trans_kg
		real*8, dimension(:), allocatable :: o2_atom_masses_kg, o_mass_kg
		real*8 :: o2_mass_kg, o_mass_sum_kg

		o2_atom_masses_kg = get_atom_masses(o2_molecule)
		o2_mass_kg = sum(o2_atom_masses_kg)

		o_mass_kg = get_atom_masses(o_atom)
		o_mass_sum_kg = sum(o_mass_kg)

		mu_trans_kg = o2_mass_kg * o_mass_sum_kg / (o2_mass_kg + o_mass_sum_kg)
	end function get_mu_trans
		
end PROGRAM RK_solution