PROGRAM RK_Solution
! Solve the differential equation dC/dt = ..  using rk4 subroutine
	implicit none
	real*8 :: t0, t_final, t, h
	real*8, dimension(:), allocatable:: C, dCdt, Cout, kdecs_per_s, C_separate, dCdt_separate, C_equilibrium_separate, part_funcs
	character(100) :: output, output1, output2, output3, output4, output5, output6, output7
	integer :: iunit, junit, kunit, nunit, munit, aunit, bunit, cunit, j, unit, n, last_bound, k
	integer :: unit_K0, unit_K1, unit_K2, unit_K3, unit_K4, unit_K5, unit_K6, unit_K7, &
	unit_K8, unit_K9, unit_K10, unit_K11, unit_K12, unit_K13, unit_K14, unit_K15, unit_K16, &
	unit_K17, unit_K18, unit_K19, unit_K20
	
! Parameters for kinetics
	character(len=100) :: directory, filename 
	character(len=100) :: directory_J0_K0, directory_J24_K0, directory_J24_K1, directory_J24_K2, &
		directory_J24_K3, directory_J24_K4, directory_J24_K5, directory_J24_K6, directory_J24_K7, directory_J24_K8, &
		directory_J24_K9, directory_J24_K10, directory_J24_K11, directory_J24_K12, directory_J24_K13, directory_J24_K14, &
		directory_J24_K15, directory_J24_K16, directory_J24_K17, directory_J24_K18, directory_J24_K19, directory_J24_K20
	
	character(len=200) :: filepath, filepath_K0, filepath_K1, filepath_K2, filepath_K3, filepath_K4, filepath_K5, filepath_K6, &
	filepath_K7, filepath_K8, filepath_K9, filepath_K10, filepath_K11, filepath_K12, filepath_K13, filepath_K14, filepath_K15, &
	filepath_K16, filepath_K17, filepath_K18, filepath_K19, filepath_K20
	
	integer :: num_states, i, ios, num_counter, num_limit, ios_K0, ios_K1, ios_K2, ios_K3, ios_K4, ios_K5, ios_K6, ios_K7, ios_K8, &
	ios_K9, ios_K10, ios_K11, ios_K12, ios_K13, ios_K14, ios_K15, ios_K16, ios_K17, ios_K18, ios_K19, ios_K20
	
	real*8, dimension(:), allocatable :: Energies, Gammas, Covalent_All, Vdw_All, Infinity_All, equilibrium_constants_m3
	real*8, dimension(:, :), allocatable :: transition_matrix
	real*8, dimension(:), allocatable :: First_term, Second_term, sum_rows, threshold_Energies_K
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
	integer :: threshold_j, Js, Ks, vib_sym_well, iteration_counter, print_freq, K_initial, K_final, band_width
	integer, dimension(:), allocatable :: num_states_K, K_value, threshold_j_values, Resonance, num_counter_K, &
	unit_K, ios_K, istart, ifinish
	
	real*8 :: start_time, end_time, end_time1, end_time2, end_time3, tmp_energies, tmp_gammas, tmp_covalent_all, &
	tmp_vdw_all, tmp_infinity_all, tmp_k_value, tmp_resonance, kt_energy_j, dE_down, kij_min, sum_rms, rms
	
	integer :: Ks_indep, cutoff_count, cutoff_count_max, K_dependent
	
	logical :: print_detail, truncation
	
!---------------------------------------------------------------------------------------------------------------------------	
	
	call cpu_time(start_time)
	
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

! Pressure related parameters
	ref_pressure_per_m3 = 6.44e+24
	pressure_ratio = M_per_m3 / ref_pressure_per_m3	
	
	temp_k = 298
	sigma0_m2 = 2300 * m_per_a0**2

! Parameters for energy transfer model	
	dE_down = -43.13
	dE(1) = -43.13 * j_per_cm
	!dE(2) = get_dE_up(dE(1), temp_k)
	
	j_per_k = get_j_per_k()
	kt_energy_j = temp_k * j_per_k;
	
! Energy Spectrum Parameters in wave numbers	
	lower_energy_lim = -10
	upper_energy_lim = 10
	upper_gamma_lim = 10
	lower_gamma_lim = 1.d0 ! 10**(-2)

! Ozone isotope Labels	
	o3_molecule = "666"
	o2_molecule = "66"
	o_atom = "6"

! Computing k0
	k0_m3_per_s = get_k0_2(o3_molecule, temp_k, sigma0_m2)

! Initial and final time and time step
	t0 = 0.0d0
! 0.01*1000e-9 / pressure_ratio = 1e-12 in case of 10000std of M_per_m3
	t_final = 135e-13 !1.35*0.01*1000e-9 / pressure_ratio
	h = 100e-19 !1*1e-9 / (pressure_ratio)
	band_width = 10
	print_freq = 10
	print_detail = .True.
	truncation = .False. !.True.
	K_dependent = 1

! Specify the directory and file name  
    directory_J0_K0 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/half_integers/J_0/K_0/symmetry_1"
	
	directory_J24_K0 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/half_integers/J_24/K_0/symmetry_1"
	
	directory_J24_K1 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/J_24/K_1/symmetry_1"
	directory_J24_K2 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/half_integers/J_24/K_2/symmetry_1"
	directory_J24_K3 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/J_24/K_3/symmetry_1"
	directory_J24_K4 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/half_integers/J_24/K_4/symmetry_1"

	directory_J24_K5 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/J_24/K_5/symmetry_1"
	directory_J24_K6 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/half_integers/J_24/K_6/symmetry_1"
	directory_J24_K7 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/J_24/K_7/symmetry_1"
	directory_J24_K8 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/half_integers/J_24/K_8/symmetry_1"
	
	directory_J24_K9 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/J_24/K_9/symmetry_1"
	directory_J24_K10 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/half_integers/J_24/K_10/symmetry_1"
	directory_J24_K11 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/J_24/K_11/symmetry_1"
	directory_J24_K12 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/half_integers/J_24/K_12/symmetry_1"

	directory_J24_K13 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/J_24/K_13/symmetry_1"
	directory_J24_K14 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/half_integers/J_24/K_14/symmetry_1"
	directory_J24_K15 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/J_24/K_15/symmetry_1"
	directory_J24_K16 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/half_integers/J_24/K_16/symmetry_1"
	
	directory_J24_K17 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/J_24/K_17/symmetry_1"
	directory_J24_K18 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/half_integers/J_24/K_18/symmetry_1"
	directory_J24_K19 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/J_24/K_19/symmetry_1"
	directory_J24_K20 = "/mmfs1/home/3436yermeka/ozone_kinetics/data/resonances/mol_666/half_integers/J_24/K_20/symmetry_1"
	
	filename = "state_properties.fwc"
	
	Js = 24
	vib_sym_well = 0
	num_states = 0
	
	K_initial = 0
	K_final = 20
	allocate(num_states_K(K_final+1))
	allocate(threshold_Energies_K(K_final+1))
	allocate(threshold_j_values(K_final+1))
	num_states_K = 0

! First loop over all K's inside one J value to count number of states	
	do Ks = K_initial, K_final
	
		if (Ks == 0) then
			directory = directory_J24_K0
		else if (Ks==1) then
			directory = directory_J24_K1
		else if (Ks==2) then
			directory = directory_J24_K2
		else if (Ks==3) then
			directory = directory_J24_K3
		else if (Ks==4) then
			directory = directory_J24_K4
		else if (Ks==5) then
			directory = directory_J24_K5
		else if (Ks==6) then
			directory = directory_J24_K6
		else if (Ks==7) then
			directory = directory_J24_K7
		else if (Ks==8) then
			directory = directory_J24_K8
		else if (Ks==9) then
			directory = directory_J24_K9
		else if (Ks==10) then
			directory = directory_J24_K10
		else if (Ks==11) then
			directory = directory_J24_K11
		else if (Ks==12) then
			directory = directory_J24_K12
		else if (Ks==13) then
			directory = directory_J24_K13
		else if (Ks==14) then
			directory = directory_J24_K14
		else if (Ks==15) then
			directory = directory_J24_K15
		else if (Ks==16) then
			directory = directory_J24_K16
		else if (Ks==17) then
			directory = directory_J24_K17
		else if (Ks==18) then
			directory = directory_J24_K18
		else if (Ks==19) then
			directory = directory_J24_K19
		else if (Ks==20) then
			directory = directory_J24_K20	
		end if
		
! Computes threshold energy
		!threshold_energy_j = get_higher_barrier_threshold()
		call get_threshold_energy_K(o3_molecule, o2_molecule, Ks, threshold_E, threshold_j)
		threshold_Energies_K(Ks+1) = threshold_E
		threshold_j_values(Ks+1) = threshold_j
		
! Construct the full file path
		filepath = trim(directory) // '/' // trim(filename)
! Open the file and skip the first line with kinetics data
		open(newunit=unit, file=filepath, status='old', action='read', iostat=ios)
		read(unit, *)
! Count the number of data points in the file
		do
			read(unit, *, iostat=ios) Energy_tmp, Total_Gamma_tmp, Covalent_All_tmp, Vdw_All_tmp, Infinity_All_tmp
			if (ios /= 0) exit
			if (Energy_tmp > threshold_Energies_K(Ks+1) + lower_energy_lim .and. Energy_tmp < threshold_Energies_K(Ks+1) + upper_energy_lim &
				.and. Total_Gamma_tmp <= upper_gamma_lim ) then
				num_states = num_states + 1
				num_states_K(Ks+1) = num_states_K(Ks+1) + 1
			end if
		end do
		write(*, *) "Number of states with K = ", Ks, "is: ", num_states_K(Ks+1)		
		close(unit)
	end do

!Print the total number of states from all data files	
		write(*, *) "Total number of states in" , K_final+1, "K blocks is: ", num_states
		
		call cpu_time(end_time1)
		write(*, *) "Time for the first reading:", end_time1-start_time, "seconds"
		
! Allocate arrays to store the filtered data
	allocate(Energies(num_states))
	allocate(Gammas(num_states))
	allocate(Covalent_All(num_states))
	allocate(Vdw_All(num_states))
	allocate(Infinity_All(num_states))
	allocate(K_value(num_states))
	allocate(Resonance(num_states))
	
	Energies = 0
	Gammas = 0
	Covalent_All = 0
	Vdw_All = 0
	Infinity_All = 0
	K_value = 0
	Resonance = 0

! Temporary variables for sorting	
	tmp_energies = 0
	tmp_gammas = 0
	tmp_covalent_all = 0
	tmp_vdw_all = 0
	tmp_infinity_all = 0
	tmp_k_value = 0
	tmp_resonance = 0

! Go through the data files again to read and store only the data we need						
		allocate(num_counter_K(K_final+1))
		num_counter = 1
		num_counter_K = 1
		allocate(unit_K(K_final+1))
		unit_K(1) = unit_K0
		unit_K(2) = unit_K1
		unit_K(3) = unit_K2
		unit_K(4) = unit_K3
		unit_K(5) = unit_K4
		
		unit_K(6) = unit_K5
		unit_K(7) = unit_K6
		unit_K(8) = unit_K7
		unit_K(9) = unit_K8
		unit_K(10) = unit_K9
		
		unit_K(11) = unit_K10
		unit_K(12) = unit_K11
		unit_K(13) = unit_K12
		unit_K(14) = unit_K13
		unit_K(15) = unit_K14
		
		unit_K(16) = unit_K15
		unit_K(17) = unit_K16
		unit_K(18) = unit_K17
		unit_K(19) = unit_K18
		unit_K(20) = unit_K19
		unit_K(21) = unit_K20
		
		allocate(ios_K(K_final+1))
		ios_K(1) = ios_K0
		ios_K(2) = ios_K1
		ios_K(3) = ios_K2
		ios_K(4) = ios_K3
		ios_K(5) = ios_K4
		
		ios_K(6) = ios_K5
		ios_K(7) = ios_K6
		ios_K(8) = ios_K7
		ios_K(9) = ios_K8
		ios_K(10) = ios_K9
		
		ios_K(11) = ios_K10
		ios_K(12) = ios_K11
		ios_K(13) = ios_K12
		ios_K(14) = ios_K13
		ios_K(15) = ios_K14
		
		ios_K(16) = ios_K15
		ios_K(17) = ios_K16
		ios_K(18) = ios_K17
		ios_K(19) = ios_K18
		ios_K(20) = ios_K19
		ios_K(21) = ios_K20
		
! Computes threshold energy
        !threshold_energy_j = get_higher_barrier_threshold()	
		filepath_K0 = trim(directory_J24_K0) // '/' // trim(filename)
		filepath_K1 = trim(directory_J24_K1) // '/' // trim(filename)
		filepath_K2 = trim(directory_J24_K2) // '/' // trim(filename)
		filepath_K3 = trim(directory_J24_K3) // '/' // trim(filename)
		filepath_K4 = trim(directory_J24_K4) // '/' // trim(filename)
		
		filepath_K5 = trim(directory_J24_K5) // '/' // trim(filename)
		filepath_K6 = trim(directory_J24_K6) // '/' // trim(filename)
		filepath_K7 = trim(directory_J24_K7) // '/' // trim(filename)
		filepath_K8 = trim(directory_J24_K8) // '/' // trim(filename)
		filepath_K9 = trim(directory_J24_K9) // '/' // trim(filename)
		
		filepath_K10 = trim(directory_J24_K10) // '/' // trim(filename)
		filepath_K11 = trim(directory_J24_K11) // '/' // trim(filename)
		filepath_K12 = trim(directory_J24_K12) // '/' // trim(filename)
		filepath_K13 = trim(directory_J24_K13) // '/' // trim(filename)
		filepath_K14 = trim(directory_J24_K14) // '/' // trim(filename)
		
		filepath_K15 = trim(directory_J24_K15) // '/' // trim(filename)
		filepath_K16 = trim(directory_J24_K16) // '/' // trim(filename)
		filepath_K17 = trim(directory_J24_K17) // '/' // trim(filename)
		filepath_K18 = trim(directory_J24_K18) // '/' // trim(filename)
		filepath_K19 = trim(directory_J24_K19) // '/' // trim(filename)
		filepath_K20 = trim(directory_J24_K20) // '/' // trim(filename)
		
		open(newunit=unit_K(1), file=filepath_K0, status='old', action='read', iostat=ios_K(1))
		read(unit_K(1), *)
		open(newunit=unit_K(2), file=filepath_K1, status='old', action='read', iostat=ios_K(2))
		read(unit_K(2), *)
		open(newunit=unit_K(3), file=filepath_K2, status='old', action='read', iostat=ios_K(3))
		read(unit_K(3), *)
		open(newunit=unit_K(4), file=filepath_K3, status='old', action='read', iostat=ios_K(4))
		read(unit_K(4), *)
		open(newunit=unit_K(5), file=filepath_K4, status='old', action='read', iostat=ios_K(5))
		read(unit_K(5), *)

		open(newunit=unit_K(6), file=filepath_K5, status='old', action='read', iostat=ios_K(6))
		read(unit_K(6), *)
		open(newunit=unit_K(7), file=filepath_K6, status='old', action='read', iostat=ios_K(7))
		read(unit_K(7), *)
		open(newunit=unit_K(8), file=filepath_K7, status='old', action='read', iostat=ios_K(8))
		read(unit_K(8), *)
		open(newunit=unit_K(9), file=filepath_K8, status='old', action='read', iostat=ios_K(9))
		read(unit_K(9), *)
		open(newunit=unit_K(10), file=filepath_K9, status='old', action='read', iostat=ios_K(10))
		read(unit_K(10), *)
		
		open(newunit=unit_K(11), file=filepath_K10, status='old', action='read', iostat=ios_K(11))
		read(unit_K(11), *)
		open(newunit=unit_K(12), file=filepath_K11, status='old', action='read', iostat=ios_K(12))
		read(unit_K(12), *)
		open(newunit=unit_K(13), file=filepath_K12, status='old', action='read', iostat=ios_K(13))
		read(unit_K(13), *)
		open(newunit=unit_K(14), file=filepath_K13, status='old', action='read', iostat=ios_K(14))
		read(unit_K(14), *)
		open(newunit=unit_K(15), file=filepath_K14, status='old', action='read', iostat=ios_K(15))
		read(unit_K(15), *)
		
		open(newunit=unit_K(16), file=filepath_K15, status='old', action='read', iostat=ios_K(16))
		read(unit_K(16), *)
		open(newunit=unit_K(17), file=filepath_K16, status='old', action='read', iostat=ios_K(17))
		read(unit_K(17), *)
		open(newunit=unit_K(18), file=filepath_K17, status='old', action='read', iostat=ios_K(18))
		read(unit_K(18), *)
		open(newunit=unit_K(19), file=filepath_K18, status='old', action='read', iostat=ios_K(19))
		read(unit_K(19), *)
		open(newunit=unit_K(20), file=filepath_K19, status='old', action='read', iostat=ios_K(20))
		read(unit_K(20), *)
		open(newunit=unit_K(21), file=filepath_K20, status='old', action='read', iostat=ios_K(21))
		read(unit_K(21), *)

! Read and store the filtered data
		do 
			do Ks = K_initial, K_final		
				do 
					if (num_counter_K(Ks+1) > num_states_K(Ks+1)) exit
					read(unit_K(Ks+1), *, iostat=ios_K(Ks+1)) Energies(num_counter), Gammas(num_counter), &
					Covalent_All(num_counter), Vdw_All(num_counter), Infinity_All(num_counter)
					if (Energies(num_counter) > threshold_Energies_K(Ks+1) + lower_energy_lim .and. &
						Energies(num_counter) < threshold_Energies_K(Ks+1) + upper_energy_lim &
								.and. Gammas(num_counter) <= upper_gamma_lim ) then
		! Neglect Gamma if state is below threshold, otherwise call it a resonance					
						if (Energies(num_counter) < threshold_Energies_K(Ks+1)) then
							Gammas(num_counter) = 0
						else 
							Resonance(num_counter) = 1
						end if
						K_value(num_counter) =  Ks! Change this in the future
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
									! Sorting K_value array
									tmp_k_value = K_value(k)
									K_value(k) = K_value(k-1)
									K_value(k-1) = tmp_k_value
									! Sorting Resonance array
									tmp_resonance = Resonance(k)
									Resonance(k) = Resonance(k-1)
									Resonance(k-1) = tmp_resonance	
								else
									exit
								end if						
							end do
						end if
						num_counter_K(Ks+1) = num_counter_K(Ks+1) + 1
						num_counter = num_counter + 1
						if (num_counter > num_states) goto 10
						exit
						end if
				end do
			end do		
	end do
		
10		close(unit_K0)
		close(unit_K1)
		close(unit_K2)
		close(unit_K3)
		close(unit_K4)
		
		close(unit_K5)
		close(unit_K6)
		close(unit_K7)
		close(unit_K8)
		close(unit_K9)
		
		close(unit_K10)
		close(unit_K11)
		close(unit_K12)
		close(unit_K13)
		close(unit_K14)

		close(unit_K15)
		close(unit_K16)
		close(unit_K17)
		close(unit_K18)
		close(unit_K19)
		close(unit_K20)
		
	write(*, *) "Counter over number of states:", num_counter-1
	call cpu_time(end_time2)
	write(*, *) "Time for the second reading and sorting:", end_time2-end_time1, "seconds"
	
! Check that sorting is correct	
	do i = 2, num_states
		if (Energies(i) < Energies(i-1)) then
			print *, 'error in sorting'
		end if
	end do

! Printing partition functions of O+O2 system in K blocks	
	allocate(part_funcs(K_final + 1))
	part_funcs = 0
	do Ks = K_initial, K_final
	
		if (K_dependent == 0) then
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
	!num_counter = 1

	do  num_counter = 1, num_states
		Ks = K_value(num_counter)
		
		if (K_dependent == 0) then
			Ks_indep = 0
		else
			Ks_indep = Ks
		end if
		
		call get_threshold_energy_K(o3_molecule, o2_molecule, Ks_indep, threshold_E, threshold_j)
		J_rot_start = threshold_j
		part_funcs_o2_per_m3 = calc_part_func_O2_per_m3_total_mol(temp_k, o2_molecule, o_atom, J_rot_start, Ks)
		
		equilibrium_constants_m3(num_counter) = calculate_formation_decay_equilibrium(Energies(num_counter), &
		temp_k, part_funcs_o2_per_m3, Js, Ks)
		print *, equilibrium_constants_m3(num_counter)
		kdecs_per_s(num_counter) = Gammas(num_counter) * j_per_cm / hbar_js
		!num_counter = num_counter + 1
		!if (num_counter > num_states) exit
	end do
	
	output = 'spectrum_info.txt'
	open(newunit=iunit, file=output, status='replace')
	do j = 1, num_states
!		if (K_value(j) == 0) then
		write(iunit, *) j, Energies(j), Gammas(j), Js, K_value(j), Resonance(j), Covalent_All(j), &
			Vdw_All(j), Infinity_All(j), equilibrium_constants_m3(j)
!		end if
	end do
	close(iunit)

		
! Sum of equilibrium constants
	K_eqs_tot = 0
	do  i = 1, num_states
		K_eqs_tot = K_eqs_tot + equilibrium_constants_m3(i)
	end do

! Compute and print transition matrix	
	allocate(transition_matrix(num_states, num_states))
    transition_matrix = calculate_transition_matrix_unitless(Energies, Covalent_All, Vdw_All, Infinity_All, dE_down)
	
	output2 = 'transition_matrix.csv'
    open(newunit=kunit, file=output2, status='replace')
	do i = 1, size(transition_matrix, 1)
	 write(kunit, '(1x,I8,a1)', advance='no') i, ','
	end do
	write (kunit, *)
	
	do i = 1, size(transition_matrix, 1)
		do j = 1, size(transition_matrix, 2)
		 if (j == 1) then
			write(kunit, '(I8,a1,1x,E19.12,a1,1x)', advance='no') i, ',', transition_matrix(i, j), ','
		 else 
			write(kunit, '(E19.12,a1,1x)', advance='no') transition_matrix(i, j), ',' 
		 end if
		!write(kunit, *) transition_matrix(i, :), ','
		end do
		write (kunit, *)
	end do
	close(kunit)

! Counting the elements which are greater than cuttoff value in the transition matrix
	output7 = 'cuttoff_count_numbers.txt'
	open(newunit=cunit, file=output7, status='replace')
	cutoff_count_max = 0
	kij_min = 0.001
	
	allocate(istart(num_states))
	allocate(ifinish(num_states))
	
	ifinish = num_states
	
	do i = 1, num_states
	  cutoff_count = 0
	  istart(i) = 0
		do j = 1, num_states
		 if (transition_matrix(j, i) .gt. kij_min) then 
		 cutoff_count = cutoff_count + 1
		 if (istart(i) .eq. 0) istart(i) = j
		 ifinish(i) = j
		 end if
		end do
		write (cunit, *) i, cutoff_count, istart(i), ifinish(i), ifinish(i) - istart(i) + 1, ifinish(i) - istart(i) + 1 - cutoff_count
		
		if (cutoff_count .gt. cutoff_count_max) cutoff_count_max = cutoff_count
	end do	
	write (cunit, *) cutoff_count_max
	close(cunit)
	
	transition_matrix = transition_matrix * k0_m3_per_s * M_per_m3

! Compute sum of rows in the transition matrix
	allocate(sum_rows(num_states))
	sum_rows = 0
	do i = 1, num_states
		do j = 1, num_states
				sum_rows(i) = sum_rows(i) + transition_matrix(i, j)
		end do
	end do
	
!  Compute Constant terms for Master-equations
	allocate(First_term(num_states))
	allocate(Second_term(num_states)) 	
	do  i = 1, num_states
		First_term(i) = kdecs_per_s(i)*equilibrium_constants_m3(i)*C_O_per_m3*C_O2_per_m3
		Second_term(i) = kdecs_per_s(i) + sum_rows(i)
	end do
	
! Setting up and printing the initial conditions for concentrations at t=0
	allocate(C(num_states))
	allocate(dCdt(num_states))
	
	t = t0
	C = 3559735.99047806 !0.d0
	
	output1 = 'propagated_concentration.txt'
	open(newunit=junit, file=output1, status='replace')
	if (print_detail .eqv. .True.) then
		write(junit, *) t, C
	end if
	
	iteration_counter = 0

! Compute Initial dCdt, Total dCdt and krec
	call derivs(t, C, dCdt)
	C_tot = 0
	dCdt_tot = 0
		do  i = 1, num_states
		   C_tot = C_tot + C(i)
		   dCdt_tot = dCdt_tot + (Covalent_All(i) + Vdw_All(i))*dCdt(i)
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
		
	output3 = 'recombination_coefficient_and_dCdt_tot.txt'
	open(newunit=nunit, file=output3, status='replace')
	write(nunit, *) t, k_rec, dCdt_tot, dCdt_separate
		
	output4 = 'total_concentrations.txt'
	open(newunit=munit, file=output4, status='replace')
	write(munit, *) t, C_tot, C_separate
	
	output5 = 'equilibrium_concentrations'
	open(newunit=aunit, file=output5, status='replace')
	write(aunit, *) h, 	 equilibrium_constants_m3 * C_O_per_m3 * C_O2_per_m3
	write(aunit, *) t_final, equilibrium_constants_m3 * C_O_per_m3 * C_O2_per_m3
	close(aunit)
	
	output6 = 'equilibrium_concentrations_in_K_blocks.txt'
	open(newunit=bunit, file=output6, status='replace')
	write(bunit, *) h, 1, C_equilibrium_separate
	write(bunit, *) t_final, 1, part_funcs
	close(bunit)
	
	call cpu_time(end_time3)
	write(*, *) "All parameters before main propagation loop:", end_time3-end_time2, "seconds"	
	
!  The main time-propogation loop
	allocate(Cout(num_states))
	do while (t<=t_final)
			  call derivs(t, C, dCdt)	  
! Solve the differential equation using rk4 subroutine and update variables
			  call rk4(C, dCdt, num_states, t, h, Cout, derivs)
			  t = t + h
			  C = Cout 	  
			  iteration_counter = iteration_counter + 1
			  
! Processing and printing of the results: Total dCdt and krec
			  if (mod(iteration_counter, print_freq) == 0) then			  
				  if (print_detail .eqv. .True.) then
					write(junit, *) t, C
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
				  
				  write(nunit, *) t, k_rec, dCdt_tot, dCdt_separate
				  write(munit, *) t, C_tot, C_separate
				  flush(nunit)
				  flush(munit)
			  end if
	
	end do
	
	close(junit)
	close(nunit)
	close(munit)
		
	call cpu_time(end_time)
	write(*, *) "CPU Time:", end_time-start_time, "seconds"

! Main part of the code ends here
!-------------------------------------------------------------------------------------------------------------
	
	contains
		 	
	subroutine derivs(t, C, dCdt)
		real*8, intent(in) :: t, C(num_states)
		real*8, intent(out) :: dCdt(num_states)
		real*8 Third_term
		
! Calculate the derivatives dC/dt
! Loop over final states, same as equation numbers:
		do i = 1, num_states
! calculate sum for individual state - Third_term		
			Third_term = 0 ! Initialize Third_term for this state	
				
! Loop over initial states, all terms of the summ:
				if (truncation .eqv. .False.) then
					do j = 1, num_states 
						Third_term = Third_term + transition_matrix(j, i) * C(j)
					end do
				else
!					istart(i) = max(1, i-band_width/2)
!					ifinish(i) = min(i+band_width/2, num_states)
					do j = istart(i), ifinish(i) 
						Third_term = Third_term + transition_matrix(j, i) * C(j)
					end do
				end if
								
			dCdt(i) = First_term(i) - Second_term(i)*C(i) + Third_term
		end do
			!dCdt = First_term - kdecs_per_s*C + MATMUL(modified_transition_matrix_per_s_transposed, C)
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
		real*8 :: ptran
		real*8, allocatable, dimension(:, :) :: matrix
		integer :: i, j
		
		allocate(matrix(num_states, num_states))
		matrix = 0
		
		do j = 1, size(matrix, 2)
		   do i = 1, j-1
			  ptran = Cov(i)*Cov(j) + Vdw(i)*Vdw(j) + Inf(i)*Inf(j)
			  dE_up = get_dE_up(dE_down, temp_k, Energies(i), Energies(j), equilibrium_constants_m3(i), equilibrium_constants_m3(j))
			  if (Energies(j) > Energies(i)) then
				matrix(i, j) = ptran * exp((Energies(i) - Energies(j)) / dE_up)
				matrix(j, i) = ptran * exp((Energies(j) - Energies(i)) / dE_down)
			  else if (Energies(j) < Energies(i)) then
				matrix(i, j) = ptran * exp((Energies(i) - Energies(j)) / dE_down)
				matrix(j, i) = ptran * exp((Energies(j) - Energies(i)) / dE_up)
			  end if 
		   end do
		end do
	
	end function calculate_transition_matrix_unitless
	   
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
		
!		call get_threshold_energy_K(o3_molecule, o2_molecule, K_value, threshold_E, threshold_j)
		j_per_k = get_j_per_k()
		kt_energy_j = temp_k * j_per_k;
		Total_Energy_j = Energy * j_per_cm
!		threshold_energy_j = threshold_E * j_per_cm
			
		Kfds_m3 = (2*J_value+1)*exp(-(Total_Energy_j) / kt_energy_j) / part_funcs_o2_per_m3
		
	end function calculate_formation_decay_equilibrium
	
	function get_dE_up(dE_down, temp_k, Energy_i, Energy_j, K_eq_i, K_eq_j) result(dE_up)
! Calculates dE_up corresponding to dE_down to satisfy the reversibility principle
		real*8, intent(in) :: dE_down, temp_k, Energy_i, Energy_j, K_eq_i, K_eq_j 
		real*8 :: j_per_k, kt_energy_j, dE_up

		j_per_k = get_j_per_k()
		kt_energy_j = temp_k * j_per_k
		!dE_up_j = dE_down_j / (dE_down_j / kt_energy_j - 1)
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
!		call get_threshold_energy_K(o3_molecule, o2_molecule, K, threshold_E, threshold_j)
!		threshold_energy_j = threshold_E * j_per_cm
		
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