program compute_intensity
    implicit none
    integer, parameter                  :: n_bins = 50000, n_el_max = 5, n_tab = 50000, n_dims = 3
    double precision, parameter         :: tab_width = 5000.0d0/n_tab, d_max = 20.0d0
    double complex, parameter           :: imag = (0.0d0, 1.0d0)
    logical, parameter                  :: use_debye_waller = .false., precise = .true., fast = .false.
    integer                             :: n_frames, i_at, j_at, i_bin, q_i, p_type, i_frame, n_atoms_cst, n_el, i_el, j_el, error
    integer, allocatable                :: typ(:, :), n_atoms(:), hist(:, :)
    double precision                    :: dr(n_dims), sinc_table(0:n_tab), f_abs2, intensity_q
    double precision, allocatable       :: xyzt(:, :, :), mean_xyz(:, :), msd(:), q_vals(:), f0(:, :), f1_f2(:, :), intensity(:),&
                                            & box(:, :)
    double complex                      :: f1, f2, f1_bare
    character(len=2)                    :: el(n_el_max)
    character(len=20)                   :: energy

    call read_movie("movie.xyz", xyzt, n_atoms, box, typ, el, n_el)
    xyzt = xyzt/10.0d0          ! A to nm
    box = box/10.0d0
    n_frames = size(xyzt, dim=3)

    call read_qvals("q_vals.dat", q_vals)

    call get_command_argument(1, energy, status=error)
    if (error.ne.0) stop "Error: Call program with x-ray photon energy in eV as argument."

!    write(*, *) "number of frames:", n_frames
!    write(*, *) "average number of atoms:", sum(n_atoms)/n_frames
!    write(*, *) "number of q values:", size(q_vals)

    allocate(f0(n_el, size(q_vals)), f1_f2(n_el, 2))
    call read_f(energy, el(1:n_el), q_vals, f0, f1_f2)

    if (use_debye_waller) then
    ! This is fast for long movies, but only an approximation, as it assumes that atoms move in an isotropic and harmonic manner.
        if (.not.all(n_atoms.eq.n_atoms(1))) stop "Error: Cannot use Debye-Waller approach, number of atoms not constant."
        write(*, *) "Warning: Debye-Waller approach assumes isotropic harmonic motion and requires box "//&
                    &"dimensions to be fixed and ordering of atoms to remain constant."
        n_atoms_cst = n_atoms(1)
        call tabulate_sinc()
        allocate(mean_xyz(n_dims, n_atoms_cst), msd(n_atoms_cst))

        mean_xyz = sum(xyzt, dim=3) / n_frames
        msd = sum(sum((xyzt - spread(mean_xyz, dim=3, ncopies=n_frames))**2, dim=1), dim=2) / (n_frames * n_dims)

        open(15, file="intensity_dbf.dat", action="write", status="replace")
        do q_i=1, size(q_vals)
            intensity_q = 0.0d0
            do i_at = 1, n_atoms_cst
                f1_bare = f1_f2(typ(i_at, 1), 1) + imag*f1_f2(typ(i_at, 1), 2) + f0(typ(i_at, 1), q_i)
                f1 = f1_bare * exp(-q_vals(q_i)**2 * msd(i_at) / 2)
                do j_at = i_at+1, n_atoms_cst
                    f2 = (f1_f2(typ(j_at, 1), 1) + imag*f1_f2(typ(j_at, 1), 2) + f0(typ(j_at, 1), q_i)) &
                            &* exp(-q_vals(q_i)**2 * msd(j_at) / 2)
                    dr = mean_xyz(:, j_at) - mean_xyz(:, i_at)
                    dr = dr - box(:, 1) * nint(dr/box(:, 1))
                    intensity_q = intensity_q + dble(f1 * conjg(f2) + f2 * conjg(f1)) * sinc_tb(q_vals(q_i)*norm2(dr))
                enddo
                intensity_q = intensity_q + dble(f1_bare * conjg(f1_bare))
            enddo
            write(15, *) q_vals(q_i), intensity_q
        enddo
        close(15)
    endif


    if(fast) then
    ! This is fast for systems with many atoms, but approximate since it relies on tabulated distance values.
        write(*, *) "Warning: Fast approach considers only distances up to d_max = ", d_max, " nm."
        allocate(hist((n_el * (n_el+1)) /2, -1:n_bins), intensity(size(q_vals)))
        call tabulate_sinc()

        intensity = 0.0d0
        do i_frame = 1, n_frames
            hist = 0
            do i_at = 1, n_atoms(i_frame)
                do j_at = i_at+1, n_atoms(i_frame)
                    dr = xyzt(:, j_at, i_frame) - xyzt(:, i_at, i_frame)
                    dr = dr - box(:, i_frame) * nint(dr/box(:, i_frame))
                    i_bin = nint(norm2(dr)/d_max*n_bins)
                    if(i_bin.le.n_bins) then
                        p_type = p_map(typ(i_at, i_frame), typ(j_at, i_frame), n_el)
                        hist(p_type, i_bin) = hist(p_type, i_bin) +1
                    endif
                enddo
                p_type = p_map(typ(i_at, i_frame), typ(i_at, i_frame), n_el)
                hist(p_type, -1) = hist(p_type, -1) + 1
            enddo

            do q_i=1, size(q_vals)
                do i_el = 1, n_el
                    f1 = f1_f2(i_el, 1) + imag*f1_f2(i_el, 2) + f0(i_el, q_i)
                    do j_el = i_el, n_el
                        f2 = f1_f2(j_el, 1) + imag*f1_f2(j_el, 2) + f0(j_el, q_i)
                        p_type = p_map(i_el, j_el, n_el)
                        f_abs2 = dble(f1 * conjg(f2) + f2 * conjg(f1))
                        intensity(q_i) = intensity(q_i) + f_abs2 * hist(p_type, -1) / 2.0d0
                        do i_bin = 1, n_bins
                            intensity(q_i) = intensity(q_i) + f_abs2 * hist(p_type, i_bin) * sinc_tb(q_vals(q_i)*i_bin*d_max/n_bins)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        intensity = intensity/n_frames

        open(16, file="intensity_fast.dat", action="write", status="replace")
        do q_i=1, size(q_vals)
            write(16, *) q_vals(q_i), intensity(q_i)
        enddo
        close(16)

    endif

    if(precise) then
    ! This is slow, but may serve as a reference. Code is also simpler and less error prone maybe.
        open(17, file="intensity_precise.dat", action="write", status="replace")
        do q_i=1, size(q_vals)
            intensity_q = 0.0d0
            do i_frame = 1, n_frames
                do i_at = 1, n_atoms(i_frame)
                    f1 = f1_f2(typ(i_at, i_frame), 1) + imag*f1_f2(typ(i_at, i_frame), 2) + f0(typ(i_at, i_frame), q_i)
                    do j_at = i_at+1, n_atoms(i_frame)
                        f2 = f1_f2(typ(j_at, i_frame), 1) + imag*f1_f2(typ(j_at, i_frame), 2) + f0(typ(j_at, i_frame), q_i)
                        dr = xyzt(:, j_at, i_frame) - xyzt(:, i_at, i_frame)
                        dr = dr - box(:, i_frame) * nint(dr/box(:, i_frame))
                        intensity_q = intensity_q + dble(f1 * conjg(f2) + f2 * conjg(f1)) * sinc(q_vals(q_i)*norm2(dr))
                    enddo
                    intensity_q = intensity_q + dble(f1 * conjg(f1))
                enddo
            enddo
            intensity_q = intensity_q/n_frames
            write(17, *) q_vals(q_i), intensity_q
        enddo
        close(17)
    endif

contains


    double precision elemental function sinc(x) result (s)
        implicit none
        double precision, intent(in)                :: x
        double precision, parameter                 :: eps = 1.0d-9

        if (abs(x).gt.eps) then
            s = sin(x)/x
        else
            s = 1.0d0
        endif
    endfunction sinc

    subroutine tabulate_sinc()
        implicit none
        integer                                     :: i
        do i=0, n_tab
            sinc_table(i) = sinc(i*tab_width)
        enddo
    endsubroutine tabulate_sinc

    double precision elemental function sinc_tb(x) result (s)
        implicit none
        double precision, intent(in)                :: x
        double precision                            :: frac, bin_f
        integer                                     :: bin

        bin_f = abs(x) / tab_width
        bin = floor(bin_f)
        frac = bin_f - bin

        if (bin+1.le.n_tab) then
            s = (1.0d0-frac)*sinc_table(bin) + frac*sinc_table(bin+1)
        else
            s = 0.0d0
        endif
    endfunction sinc_tb


    integer elemental function p_map(i, j, n) result (k)
        implicit none
        integer, intent(in)                         :: i, j, n
        integer                                     :: i1, j1
        i1 = max(i, j)
        j1 = min(i, j)
        k = i1 + (j1 -1) * n - (j1 * (j1-1))/2
    endfunction p_map


    subroutine read_movie(filename, xyzt, n_atoms, box, typ, el, n_el)
        implicit none
        character(len=*), intent(in)                :: filename
        double precision, allocatable, intent(out)  :: xyzt(:, :, :), box(:, :)
        integer, allocatable, intent(out)           :: n_atoms(:)
        integer, allocatable, intent(out)           :: typ(:, :)
        character(len=2), intent(out)               :: el(n_el_max)
        integer, intent(out)                        :: n_el
        integer                                     :: n_at, n_frames, error, i_at, i_frame, n_atoms_max, i_el
        character(len=2), allocatable               :: typ_c(:, :)
        logical                                     :: found

        open(20, file=trim(filename), action='read', iostat=error)
        if (error.ne.0) stop "Error in function read_movie: "//filename//" missing"

        error = 0
        n_frames = 0
        n_atoms_max = 0
        do while (.true.)
            read(20, *, iostat=error) n_at
            if (error.ne.0) exit
            if (n_at.gt.n_atoms_max) n_atoms_max = n_at
            read(20, *)
            do i_at = 1, n_at
                read(20, *)
            enddo
            n_frames = n_frames + 1
        enddo
        rewind(20)

        allocate(xyzt(n_dims, n_atoms_max, n_frames), typ_c(n_atoms_max, n_frames), n_atoms(n_frames), box(n_dims, n_frames))
        allocate(typ(n_atoms_max, n_frames))

        do i_frame = 1, n_frames
            read(20, *) n_atoms(i_frame)
            ! Read box dimensions from the xyz file if they are not infinite
            ! read(20, *, iostat=error) box(1, i_frame), box(2, i_frame), box(3, i_frame)
            read(20, *)
            box(:, i_frame) = huge(box)
            do i_at = 1, n_atoms(i_frame)
                read(20, *) typ_c(i_at, i_frame), xyzt(:, i_at, i_frame)
            enddo
        enddo
        close(20)

        ! Map the character atom types to integer element indices per frame
        el = "  "
        n_el = 0
        do i_frame = 1, n_frames
            do i_at = 1, n_atoms(i_frame)
                found = .false.
                do i_el = 1, n_el
                    if (typ_c(i_at, i_frame).eq.el(i_el)) then
                        typ(i_at, i_frame) = i_el
                        found = .true.
                        exit
                    endif
                enddo
                if (found.eqv..false.) then
                    n_el = n_el + 1
                    if (n_el.gt.n_el_max) stop "Error: increase n_el_max"
                    typ(i_at, i_frame) = n_el
                    el(n_el) = typ_c(i_at, i_frame)
                endif
            enddo
        enddo

    endsubroutine read_movie


    subroutine read_qvals(filename, q_vals)
        implicit none
        integer                                     :: i_line, n_lines, error
        character(len=*), intent(in)                :: filename
        double precision, allocatable, intent(out)  :: q_vals(:)

        open(18, file=filename, action='read', iostat=error)
        if (error.ne.0) stop "Error in function read_qvals: "//filename//" missing"

        error = 0
        n_lines = 0
        do while (error.eq.0)
            read(18, *, iostat=error)
            n_lines = n_lines + 1
        enddo
        rewind(18)
        n_lines = n_lines-1

        allocate(q_vals(n_lines))

        do i_line = 1, n_lines
            read(18, *) q_vals(i_line)
        enddo
        close(18)

    endsubroutine read_qvals


    subroutine read_f(energy, el, q_vals, f0, f1_f2)
        implicit none
        double precision, intent(in)                :: q_vals(:)
        character(len=2), intent(in)                :: el(:)
        character(len=20), intent(in)               :: energy
        double precision, intent(out)               :: f0(:, :), f1_f2(:, :)
        integer                                     :: i_line, i_el, error
        character(len=500)                          :: command_string
        character(len=64)                           :: f_file

        ! Unique per-process output name so concurrent runs do not clobber each other.
        write(f_file, '(A, I0, A)') "f_vals_", getpid(), ".dat"

        command_string = "python3 write_f.py "//trim(energy)//" "//trim(f_file)
        do i_el = 1, size(el)
            command_string = trim(command_string)//" "//trim(el(i_el))
        enddo
        call execute_command_line(trim(command_string), exitstat=error)
        if (error.ne.0) stop "Error in function read_f: python script error"

        open(19, file=trim(f_file), action="read")
        do i_line = 1, size(q_vals)
            read(19, *) f0(1:size(el), i_line)
        enddo
        do i_el=1, size(el)
            read(19, *) f1_f2(i_el, :)
        enddo
        close(19, status="delete")

    endsubroutine read_f


endprogram

