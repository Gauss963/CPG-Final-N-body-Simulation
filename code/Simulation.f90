program nbody_sim
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: exp_val, num_particles, steps, i, j, step
    real(dp) :: L, dt, T, G

    real(dp), allocatable :: mass(:), radius(:)
    real(dp), allocatable :: pos(:,:), vel(:,:), acc(:,:)
    real(dp), allocatable :: total_kinetic_energy(:), total_potential_energy(:), total_angular_momentum(:)
    integer :: unit_particle_data, unit_energy_data
    character(len=100) :: pc_filename, eg_filename

    call read_config("config.ini", exp_val, L, dt, T, G)
    num_particles = 2**exp_val
    steps = int(T/dt)

    allocate(mass(num_particles), radius(num_particles), pos(3,num_particles), vel(3,num_particles), acc(3,num_particles))
    allocate(total_kinetic_energy(steps), total_potential_energy(steps), total_angular_momentum(steps))

    call random_seed()
    call init_particles(num_particles, L, mass, radius, pos, vel)

    ! write(pc_filename,'("../data/FF_particle_data_N_",I0,".bin")') num_particles
    ! open(newunit=unit_particle_data, file=pc_filename, form='unformatted', access='stream', status='replace')

    write(eg_filename,'("../data/FF_energy_data_N_",I0,".bin")') num_particles
    open(newunit=unit_energy_data, file=eg_filename, form='unformatted', access='stream', status='replace')

    do step = 1, steps
        write(*,*) "Step ", step, "/", steps, " Progress: ", (step * 100.0 / steps)

        call compute_gravity_forces(num_particles, pos, mass, G, acc)
        call handle_collisions(num_particles, pos, vel, mass, radius)

        do i = 1, num_particles
            vel(:,i) = vel(:,i) + acc(:,i)*dt
            pos(:,i) = pos(:,i) + vel(:,i)*dt
        enddo

        total_kinetic_energy(step) = compute_kinetic_energy(num_particles, mass, vel)
        total_potential_energy(step) = compute_potential_energy(num_particles, pos, mass, G)
        total_angular_momentum(step) = compute_angular_momentum(num_particles, pos, vel, mass)

        ! do i = 1, num_particles
        !     write(unit_particle_data) pos(1,i), pos(2,i), pos(3,i), vel(1,i), vel(2,i), vel(3,i), mass(i), radius(i)
        ! enddo
    enddo

    do i = 1, steps
        write(unit_energy_data) total_kinetic_energy(i)
    enddo
    do i = 1, steps
        write(unit_energy_data) total_potential_energy(i)
    enddo
    do i = 1, steps
        write(unit_energy_data) total_angular_momentum(i)
    enddo

    ! close(unit_particle_data)
    close(unit_energy_data)

    write(*,*) "Simulation finished."

contains

    subroutine read_config(filename, exp_val, L, dt, T, G)
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(out) :: exp_val
        real(dp), intent(out) :: L, dt, T, G

        integer :: unit, ios, eqpos
        character(len=256) :: line, section, key, val
        logical :: in_simulation_section

        in_simulation_section = .false.
        exp_val = -1; L=0.0_dp; dt=0.0_dp; T=0.0_dp; G=0.0_dp

        open(newunit=unit, file=filename, form='formatted', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Could not open config file:", filename
            stop 1
        endif

        do
            read(unit,'(A)', iostat=ios) line
            if (ios /= 0) exit

            line = trim(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == ";" .or. line(1:1) == "#") cycle

            if (line(1:1) == "[" .and. line(len_trim(line):len_trim(line)) == "]") then
                section = line(2:len_trim(line)-1)
                if (section == "Simulation") then
                    in_simulation_section = .true.
                else
                    in_simulation_section = .false.
                endif
                cycle
            endif

            if (in_simulation_section) then
                eqpos = index(line, "=")
                if (eqpos > 0) then
                    key = adjustl(trim(line(1:eqpos-1)))
                    val = adjustl(trim(line(eqpos+1:len_trim(line))))

                    select case (key)
                    case("n")
                        read(val,*) exp_val
                    case("L")
                        read(val,*) L
                    case("dt")
                        read(val,*) dt
                    case("T")
                        read(val,*) T
                    case("G")
                        read(val,*) G
                    end select
                endif
            endif
        enddo

        close(unit)

        if (exp_val<0) then
            print *, "Missing n in config"
            stop 1
        endif
        if (L==0.0_dp) then
            print *, "Missing L in config"
            stop 1
        endif
        if (dt==0.0_dp) then
            print *, "Missing dt in config"
            stop 1
        endif
        if (T==0.0_dp) then
            print *, "Missing T in config"
            stop 1
        endif
        if (G==0.0_dp) then
            print *, "Missing G in config"
            stop 1
        endif

    end subroutine read_config

    subroutine init_particles(N, L, mass, radius, pos, vel)
        implicit none
        integer, intent(in) :: N
        real(dp), intent(in) :: L
        real(dp), intent(out) :: mass(N), radius(N)
        real(dp), intent(out) :: pos(3,N), vel(3,N)
        integer :: i
        real(dp) :: r1

        do i = 1, N
            call random_number(r1); mass(i) = 1.0_dp + r1*(2.0_dp-1.0_dp)
            call random_number(r1); pos(1,i) = -L/2.0_dp + r1*L
            call random_number(r1); pos(2,i) = -L/2.0_dp + r1*L
            call random_number(r1); pos(3,i) = -L/2.0_dp + r1*L

            call random_number(r1); vel(1,i) = -1.0_dp + r1*2.0_dp
            call random_number(r1); vel(2,i) = -1.0_dp + r1*2.0_dp
            call random_number(r1); vel(3,i) = -1.0_dp + r1*2.0_dp

            call random_number(r1); radius(i) = 0.001_dp + r1*(0.05_dp-0.001_dp)
        enddo
    end subroutine init_particles

    subroutine compute_gravity_forces(N, pos, mass, G, acc)
        implicit none
        integer, intent(in) :: N
        real(dp), intent(in) :: pos(3,N), mass(N), G
        real(dp), intent(out) :: acc(3,N)
        integer :: i, j
        real(dp), allocatable :: forces(:,:)
        real(dp) :: rx, ry, rz, r2, r, Fmag, invr

        allocate(forces(3,N))
        forces = 0.0_dp

        do i = 1, N
            do j = i+1, N
                rx = pos(1,j)-pos(1,i)
                ry = pos(2,j)-pos(2,i)
                rz = pos(3,j)-pos(3,i)
                r2 = rx*rx + ry*ry + rz*rz
                if (r2 > 0.0_dp) then
                    r = sqrt(r2)
                    Fmag = G*mass(i)*mass(j)/r2
                    invr = 1.0_dp/r
                    forces(1,i) = forces(1,i) + Fmag*rx*invr
                    forces(2,i) = forces(2,i) + Fmag*ry*invr
                    forces(3,i) = forces(3,i) + Fmag*rz*invr
                    forces(1,j) = forces(1,j) - Fmag*rx*invr
                    forces(2,j) = forces(2,j) - Fmag*ry*invr
                    forces(3,j) = forces(3,j) - Fmag*rz*invr
                endif
            enddo
        enddo

        do i = 1, N
            acc(:,i) = forces(:,i)/mass(i)
        enddo

        deallocate(forces)
    end subroutine compute_gravity_forces

    subroutine handle_collisions(N, pos, vel, mass, radius)
        implicit none
        integer, intent(in) :: N
        real(dp), intent(inout) :: pos(3,N), vel(3,N)
        real(dp), intent(in) :: mass(N), radius(N)
        integer :: i, j
        real(dp) :: rx, ry, rz, dist, rad_sum
        real(dp) :: vxr, vyr, vzr, vdotn, nx, ny, nz
        real(dp) :: m1, m2, Jimp

        do i = 1, N
            do j = i+1, N
                rx = pos(1,j)-pos(1,i)
                ry = pos(2,j)-pos(2,i)
                rz = pos(3,j)-pos(3,i)
                dist = sqrt(rx*rx+ry*ry+rz*rz)
                rad_sum = radius(i)+radius(j)
                if (dist <= rad_sum) then
                    if (dist == 0.0_dp) cycle
                    nx = rx/dist
                    ny = ry/dist
                    nz = rz/dist

                    vxr = vel(1,i)-vel(1,j)
                    vyr = vel(2,i)-vel(2,j)
                    vzr = vel(3,i)-vel(3,j)

                    vdotn = vxr*nx + vyr*ny + vzr*nz

                    if (vdotn > 0.0_dp) cycle

                    m1 = mass(i)
                    m2 = mass(j)
                    Jimp = (2.0_dp * vdotn) / (1.0_dp/m1 + 1.0_dp/m2)

                    vel(1,i) = vel(1,i) - (Jimp*nx)/m1
                    vel(2,i) = vel(2,i) - (Jimp*ny)/m1
                    vel(3,i) = vel(3,i) - (Jimp*nz)/m1

                    vel(1,j) = vel(1,j) + (Jimp*nx)/m2
                    vel(2,j) = vel(2,j) + (Jimp*ny)/m2
                    vel(3,j) = vel(3,j) + (Jimp*nz)/m2
                endif
            enddo
        enddo

    end subroutine handle_collisions

    function compute_kinetic_energy(N, mass, vel) result(KE)
        implicit none
        integer, intent(in) :: N
        real(dp), intent(in) :: mass(N), vel(3,N)
        real(dp) :: KE, v2
        integer :: i
        KE = 0.0_dp
        do i = 1, N
            v2 = vel(1,i)*vel(1,i)+vel(2,i)*vel(2,i)+vel(3,i)*vel(3,i)
            KE = KE + 0.5_dp*mass(i)*v2
        enddo
    end function compute_kinetic_energy

    function compute_potential_energy(N, pos, mass, G) result(PE)
        implicit none
        integer, intent(in) :: N
        real(dp), intent(in) :: pos(3,N), mass(N), G
        real(dp) :: PE
        integer :: i,j
        real(dp) :: rx, ry, rz, r2, r
        PE = 0.0_dp
        do i = 1, N
            do j = i+1, N
                rx = pos(1,j)-pos(1,i)
                ry = pos(2,j)-pos(2,i)
                rz = pos(3,j)-pos(3,i)
                r2 = rx*rx+ry*ry+rz*rz
                if (r2 > 0.0_dp) then
                    r = sqrt(r2)
                    PE = PE - G*mass(i)*mass(j)/r
                endif
            enddo
        enddo
    end function compute_potential_energy

    function compute_angular_momentum(N, pos, vel, mass) result(L_mag)
        implicit none
        integer, intent(in) :: N
        real(dp), intent(in) :: pos(3,N), vel(3,N), mass(N)
        real(dp) :: L_mag
        real(dp) :: Lx, Ly, Lz
        real(dp) :: mvx, mvy, mvz
        integer :: i
        Lx = 0.0_dp; Ly=0.0_dp; Lz=0.0_dp
        do i = 1, N
            mvx = mass(i)*vel(1,i)
            mvy = mass(i)*vel(2,i)
            mvz = mass(i)*vel(3,i)
            Lx = Lx + (pos(2,i)*mvz - pos(3,i)*mvy)
            Ly = Ly + (pos(3,i)*mvx - pos(1,i)*mvz)
            Lz = Lz + (pos(1,i)*mvy - pos(2,i)*mvx)
        enddo
        L_mag = sqrt(Lx*Lx + Ly*Ly + Lz*Lz)
    end function compute_angular_momentum

end program nbody_sim