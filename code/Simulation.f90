program Simulation
    use iso_fortran_env, only: real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    integer, parameter :: dp = real64
    integer, parameter :: leaf_capacity = 4, max_depth = 48

    type :: tree_node
        integer :: first, last, children(8)
        real(dp) :: lower(3), upper(3), center_of_mass(3)
        real(dp) :: mass, max_radius, width_squared
    end type tree_node

    integer :: exp_val, num_particles, steps, i, step
    real(dp) :: L, dt, T, G, step_count, initial_energy
    logical :: use_tree = .false.
    real(dp) :: theta = 0.5_dp
    character(len=256) :: initial_conditions = ''

    real(dp), allocatable :: mass(:), radius(:), pos(:,:), vel(:,:), acc(:,:)
    real(dp), allocatable :: total_kinetic_energy(:), total_potential_energy(:), total_angular_momentum(:)
    integer :: unit_energy_data, ios, command_status, seed_size
    character(len=100) :: eg_filename
    integer, allocatable :: seed(:)

    type(tree_node), allocatable :: nodes(:)
    integer, allocatable :: order(:), scratch(:), rank(:), contacts(:)
    integer :: node_count

    call read_config("config.ini")
    if (exp_val < 0 .or. exp_val >= digits(num_particles)) error stop 'n is outside the range for N = 2^n'
    if (.not. ieee_is_finite(L) .or. L <= 0.0_dp .or. .not. ieee_is_finite(G)) then
        error stop 'L must be positive and G must be finite'
    end if
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp .or. .not. ieee_is_finite(T) .or. T < 0.0_dp) then
        error stop 'Invalid dt or T'
    end if
    step_count = T / dt
    if (step_count > real(huge(steps), dp)) error stop 'Too many time steps'
    num_particles = 2**exp_val
    steps = int(step_count)

    allocate(mass(num_particles), radius(num_particles), pos(3,num_particles), &
             vel(3,num_particles), acc(3,num_particles), stat=ios)
    if (ios /= 0) error stop 'Not enough memory for the particles'
    allocate(total_kinetic_energy(steps), total_potential_energy(steps), total_angular_momentum(steps), stat=ios)
    if (ios /= 0) error stop 'Not enough memory for the time series'

    if (len_trim(initial_conditions) == 0) then
        call random_seed(size=seed_size)
        allocate(seed(seed_size))
        seed = 42
        call random_seed(put=seed)
        call init_particles(num_particles, L, mass, radius, pos, vel)
    else
        call read_particles(initial_conditions)
    end if

    if (use_tree) then
        write(*,'(A,F0.4,A)') 'Gravity: barnes_hut (theta=', theta, ')'
        allocate(nodes(min(num_particles, 1024)), order(num_particles), scratch(num_particles), &
                 rank(num_particles), contacts(num_particles), stat=ios)
        if (ios /= 0) error stop 'Not enough memory for the octree'
    else
        write(*,'(A)') 'Gravity: direct'
    end if

    ! The tree and accelerations from the previous evaluation belong to the current positions.
    if (steps > 0) call evaluate_gravity(initial_energy)
    do step = 1, steps
        write(*,*) 'Step ', step, '/', steps, ' Progress: ', (step * 100.0_dp / steps)

        call handle_collisions()
        do i = 1, num_particles
            vel(:,i) = vel(:,i) + acc(:,i) * dt
            pos(:,i) = pos(:,i) + vel(:,i) * dt
        end do

        total_kinetic_energy(step) = compute_kinetic_energy(num_particles, mass, vel)
        call evaluate_gravity(total_potential_energy(step))
        total_angular_momentum(step) = compute_angular_momentum(num_particles, pos, vel, mass)
    end do

    call execute_command_line('mkdir -p ../data', exitstat=ios, cmdstat=command_status)
    if (command_status /= 0) error stop 'Could not create ../data'
    if (ios /= 0) error stop 'Could not create ../data'
    write(eg_filename,'("../data/FF_energy_data_N_",I0,".bin")') num_particles
    open(newunit=unit_energy_data, file=eg_filename, form='unformatted', access='stream', status='replace', iostat=ios)
    if (ios /= 0) error stop 'Could not open energy output'
    write(unit_energy_data, iostat=ios) total_kinetic_energy, total_potential_energy, total_angular_momentum
    if (ios /= 0) error stop 'Could not write energy output'
    close(unit_energy_data, iostat=ios)
    if (ios /= 0) error stop 'Could not close energy output'

    write(*,*) "Simulation finished."

contains

    subroutine read_config(filename)
        implicit none
        character(len=*), intent(in) :: filename

        integer :: unit, ios, eqpos
        character(len=256) :: line, section, key, val
        logical :: in_simulation_section

        in_simulation_section = .false.
        exp_val = -1
        L = ieee_value(0.0_dp, ieee_quiet_nan)
        dt = L
        T = L
        G = L

        open(newunit=unit, file=filename, form='formatted', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Could not open config file:", filename
            stop 1
        endif

        do
            read(unit,'(A)', iostat=ios) line
            if (is_iostat_end(ios)) exit
            if (ios /= 0) error stop 'Could not read config.ini'

            line = adjustl(line)
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
                        read(val,*,iostat=ios) exp_val
                    case("L")
                        read(val,*,iostat=ios) L
                    case("dt")
                        read(val,*,iostat=ios) dt
                    case("T")
                        read(val,*,iostat=ios) T
                    case("G")
                        read(val,*,iostat=ios) G
                    case("gravity")
                        select case (val)
                        case ("direct")
                            use_tree = .false.
                        case ("barnes_hut")
                            use_tree = .true.
                        case default
                            error stop "gravity must be direct or barnes_hut"
                        end select
                    case("theta")
                        read(val,*,iostat=ios) theta
                    case("initial_conditions")
                        initial_conditions = val
                    end select
                    if (ios /= 0) error stop "Invalid numeric value in config.ini"
                endif
            endif
        enddo

        close(unit)

        if (.not. ieee_is_finite(theta) .or. theta < 0.0_dp) then
            error stop "theta must be a finite, non-negative number"
        end if
    end subroutine read_config

    subroutine init_particles(N, L, mass, radius, pos, vel)
        implicit none
        integer, intent(in) :: N
        real(dp), intent(in) :: L
        real(dp), intent(out) :: mass(N), radius(N)
        real(dp), intent(out) :: pos(3,N), vel(3,N)
        integer :: i
        real(dp) :: r1

        r1 = 42.0

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

    subroutine read_particles(filename)
        character(len=*), intent(in) :: filename
        integer :: unit, status, body
        integer(int64) :: bytes
        real(dp) :: record(8)

        open(newunit=unit, file=filename, access='stream', form='unformatted', status='old', action='read', iostat=status)
        if (status /= 0) error stop 'Could not open initial_conditions'
        inquire(unit=unit, size=bytes, iostat=status)
        if (status /= 0 .or. bytes /= 64_int64 * num_particles) error stop 'Expected N records of 8 doubles'
        do body = 1, num_particles
            read(unit, iostat=status) record
            if (status /= 0) error stop 'Could not read initial_conditions'
            if (.not. all(ieee_is_finite(record)) .or. record(7) <= 0.0_dp .or. record(8) < 0.0_dp) then
                error stop 'Invalid particle in initial_conditions'
            end if
            pos(:,body) = record(1:3)
            vel(:,body) = record(4:6)
            mass(body) = record(7)
            radius(body) = record(8)
        end do
        close(unit)
    end subroutine read_particles

    subroutine evaluate_gravity(potential)
        real(dp), intent(out) :: potential
        real(dp) :: particle_potential, theta_squared
        integer :: body, root

        if (.not. use_tree) then
            call direct_gravity(potential)
            return
        end if

        node_count = 0
        do body = 1, num_particles
            order(body) = body
        end do
        call build_node(1, num_particles + 1, 0, root)
        do body = 1, num_particles
            rank(order(body)) = body
        end do

        theta_squared = theta * theta
        potential = 0.0_dp
        do body = 1, num_particles
            acc(:,body) = 0.0_dp
            particle_potential = 0.0_dp
            call accumulate_gravity(root, rank(body), pos(:,body), theta_squared, acc(:,body), particle_potential)
            potential = potential + mass(body) * particle_potential
        end do
        potential = 0.5_dp * potential
    end subroutine evaluate_gravity

    subroutine direct_gravity(potential)
        real(dp), intent(out) :: potential
        real(dp) :: delta(3), direction(3), r2, inv_r, scale
        integer :: a, b

        potential = 0.0_dp
        acc = 0.0_dp
        do a = 1, num_particles
            do b = a + 1, num_particles
                delta = pos(:,b) - pos(:,a)
                r2 = dot_product(delta, delta)
                if (r2 == 0.0_dp) cycle
                inv_r = 1.0_dp / sqrt(r2)
                scale = G / r2
                direction = inv_r * delta
                acc(:,a) = acc(:,a) + (scale * mass(b)) * direction
                acc(:,b) = acc(:,b) - (scale * mass(a)) * direction
                potential = potential - G * mass(a) * mass(b) * inv_r
            end do
        end do
    end subroutine direct_gravity

    integer function octant(position, midpoint) result(oct)
        real(dp), intent(in) :: position(3), midpoint(3)
        oct = 1
        if (position(1) >= midpoint(1)) oct = oct + 1
        if (position(2) >= midpoint(2)) oct = oct + 2
        if (position(3) >= midpoint(3)) oct = oct + 4
    end function octant

    recursive subroutine build_node(first, last, depth, id)
        integer, intent(in) :: first, last, depth
        integer, intent(out) :: id
        type(tree_node) :: node
        type(tree_node), allocatable :: grown(:)
        real(dp) :: extent(3), midpoint(3), width
        integer :: p, body, offsets(9), next(9), oct, child, slot, status

        node_count = node_count + 1
        id = node_count
        if (id > size(nodes)) then
            allocate(grown(min(2_int64 * size(nodes), int(huge(id), int64))), stat=status)
            if (status /= 0) error stop 'Not enough memory for the octree'
            grown(1:size(nodes)) = nodes
            call move_alloc(grown, nodes)
        end if

        node%first = first
        node%last = last
        node%children = 0
        node%lower = pos(:,order(first))
        node%upper = node%lower
        node%mass = 0.0_dp
        node%center_of_mass = 0.0_dp
        node%max_radius = 0.0_dp
        do p = first, last - 1
            body = order(p)
            node%lower = min(node%lower, pos(:,body))
            node%upper = max(node%upper, pos(:,body))
            node%mass = node%mass + mass(body)
            node%center_of_mass = node%center_of_mass + mass(body) * pos(:,body)
            node%max_radius = max(node%max_radius, radius(body))
        end do
        extent = node%upper - node%lower
        width = maxval(extent)
        node%width_squared = width * width
        node%center_of_mass = node%center_of_mass / node%mass
        nodes(id) = node
        if (last - first <= leaf_capacity .or. width == 0.0_dp .or. depth == max_depth) return

        ! Partition indices by the octants of the tight bounding box.
        midpoint = node%lower + 0.5_dp * extent
        where (midpoint == node%lower) midpoint = node%upper
        offsets = 0
        offsets(1) = first
        do p = first, last - 1
            oct = octant(pos(:,order(p)), midpoint)
            offsets(oct + 1) = offsets(oct + 1) + 1
        end do
        do oct = 1, 8
            offsets(oct + 1) = offsets(oct + 1) + offsets(oct)
        end do
        next = offsets
        do p = first, last - 1
            body = order(p)
            oct = octant(pos(:,body), midpoint)
            scratch(next(oct)) = body
            next(oct) = next(oct) + 1
        end do
        order(first:last-1) = scratch(first:last-1)

        slot = 0
        do oct = 1, 8
            if (offsets(oct) == offsets(oct + 1)) cycle
            call build_node(offsets(oct), offsets(oct + 1), depth + 1, child)
            slot = slot + 1
            nodes(id)%children(slot) = child
        end do
    end subroutine build_node

    subroutine add_gravity(delta, r2, gm, acceleration, potential)
        real(dp), intent(in) :: delta(3), r2, gm
        real(dp), intent(inout) :: acceleration(3), potential
        real(dp) :: inv_r, scale

        if (r2 == 0.0_dp) return
        inv_r = 1.0_dp / sqrt(r2)
        scale = gm / r2
        acceleration = acceleration + scale * (inv_r * delta)
        potential = potential - gm * inv_r
    end subroutine add_gravity

    recursive subroutine accumulate_gravity(id, target, position, theta_squared, acceleration, potential)
        integer, intent(in) :: id, target
        real(dp), intent(in) :: position(3), theta_squared
        real(dp), intent(inout) :: acceleration(3), potential
        integer :: p, source, slot, child
        real(dp) :: delta(3), r2

        if (nodes(id)%children(1) == 0) then
            do p = nodes(id)%first, nodes(id)%last - 1
                if (p == target) cycle
                source = order(p)
                delta = pos(:,source) - position
                call add_gravity(delta, dot_product(delta, delta), G * mass(source), acceleration, potential)
            end do
            return
        end if

        ! A cell containing the target must open to avoid self-gravity.
        if (target < nodes(id)%first .or. target >= nodes(id)%last) then
            delta = nodes(id)%center_of_mass - position
            r2 = dot_product(delta, delta)
            if (nodes(id)%width_squared < theta_squared * r2) then
                call add_gravity(delta, r2, G * nodes(id)%mass, acceleration, potential)
                return
            end if
        end if
        do slot = 1, 8
            child = nodes(id)%children(slot)
            if (child == 0) exit
            call accumulate_gravity(child, target, position, theta_squared, acceleration, potential)
        end do
    end subroutine accumulate_gravity

    subroutine handle_collisions()
        integer :: a, b, count, p

        do a = 1, num_particles - 1
            if (use_tree) then
                count = 0
                call find_contacts(1, a, count)
                ! Apply the impulses in the same pair order as the direct solver.
                call sort_contacts(1, count)
                do p = 1, count
                    call collide(a, contacts(p))
                end do
            else
                do b = a + 1, num_particles
                    call collide(a, b)
                end do
            end if
        end do
    end subroutine handle_collisions

    recursive subroutine find_contacts(id, a, count)
        integer, intent(in) :: id, a
        integer, intent(inout) :: count
        integer :: p, b, slot, child
        real(dp) :: offset(3)

        offset = max(nodes(id)%lower - pos(:,a), pos(:,a) - nodes(id)%upper, 0.0_dp)
        if (norm(offset) > radius(a) + nodes(id)%max_radius) return
        if (nodes(id)%children(1) == 0) then
            do p = nodes(id)%first, nodes(id)%last - 1
                b = order(p)
                if (b <= a) cycle
                if (norm(pos(:,a) - pos(:,b)) <= radius(a) + radius(b)) then
                    count = count + 1
                    contacts(count) = b
                end if
            end do
            return
        end if
        do slot = 1, 8
            child = nodes(id)%children(slot)
            if (child == 0) exit
            call find_contacts(child, a, count)
        end do
    end subroutine find_contacts

    recursive subroutine sort_contacts(first, last)
        integer, intent(in) :: first, last
        integer :: lo, hi, left, right, pivot, saved

        lo = first
        hi = last
        do while (lo < hi)
            left = lo
            right = hi
            pivot = contacts(lo + (hi - lo) / 2)
            do
                do while (contacts(left) < pivot)
                    left = left + 1
                end do
                do while (contacts(right) > pivot)
                    right = right - 1
                end do
                if (left > right) exit
                saved = contacts(left)
                contacts(left) = contacts(right)
                contacts(right) = saved
                left = left + 1
                right = right - 1
                if (left > right) exit
            end do
            ! Recurse on the smaller partition to bound stack usage.
            if (right - lo < hi - left) then
                call sort_contacts(lo, right)
                lo = left
            else
                call sort_contacts(left, hi)
                hi = right
            end if
        end do
    end subroutine sort_contacts

    real(dp) function norm(vector) result(length)
        real(dp), intent(in) :: vector(3)
        length = sqrt(dot_product(vector, vector))
    end function norm

    subroutine collide(a, b)
        integer, intent(in) :: a, b
        real(dp) :: normal(3), distance, along_normal, impulse(3)

        normal = pos(:,a) - pos(:,b)
        distance = norm(normal)
        if (distance == 0.0_dp .or. distance > radius(a) + radius(b)) return
        normal = (1.0_dp / distance) * normal
        along_normal = dot_product(vel(:,a) - vel(:,b), normal)
        if (along_normal > 0.0_dp) return

        impulse = ((2.0_dp * along_normal) / (1.0_dp / mass(a) + 1.0_dp / mass(b))) * normal
        vel(:,a) = vel(:,a) - impulse / mass(a)
        vel(:,b) = vel(:,b) + impulse / mass(b)
    end subroutine collide
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

end program Simulation
