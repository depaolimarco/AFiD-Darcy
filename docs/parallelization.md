# Parallelization

## Pencil decomposition

### Transforms

## Halo cell update

The interface [`update_halo`](..src/2DECOMP/decomp_2d.F90) which calls subroutine [`update_halo_real`](..src/2DECOMP/halo.F90) is called on the array to update data in the halo cells at the boundaries of the pencils. 