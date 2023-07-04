module kpath

  use utility
  use data_structures

  implicit none

  type, extends(local_k_data) :: k_path !TODO: MAKE CONSTRUCTOR.
    real(kind=dp), allocatable :: vectors(:, :) !1st index is the index of the vector in the path. 2nd index corresponds to the component of the vector in the path, so its size is vectors(size(number_of_vecotrs), 3).
    integer, allocatable       :: number_of_vectors(:) !Its size is the number of vectors in the BZ, vector(i) contains the number of k-points between vector i and vector i+1.
  end type k_path

end module kpath