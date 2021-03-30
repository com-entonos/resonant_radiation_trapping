c   rad.inc

       integer*4 ZiD                     !max spatial cells
       parameter(ZiD=55)
       integer*4 kid                     !max frequency cells
       parameter(kid=105)

c  physical constants
       real*4 amu,eCharge,m,Kb
       real*8 pi
       parameter(Kb=1.380662e-23)
       parameter(amu=1.6605655e-27)
       parameter(eCharge=1.6021892e-19)
       parameter(m=9.109534e-31)
       parameter(pi=3.14159265359d0)

       real nra(0:zid+1,kid,1:2)              !number of excited atoms (z,f,curr)
       real prra(0:zid,kid)              !production (#/s) of excited atoms (z,f)
       real pa(0:zid+1,0:zid,kid)       !prob of going to final from initial: Q
       real pac(kid,kid)                     !freq redistr. prob: J
       real*8 kv(kid),fv(kid),df(2,kid)   !absorption, freq, freq bounds
       real*4 zmax,dz
       real*8 zmaxdp,dzdp,k0,tvac,av,pcoll
       real*8 int_ls(kid), nesc(0:1,kid)
       integer*4 ci,izmax,ikmax,ikbmax
       character*80 prodfile,filedat
       character*40 pathname
       integer*4 ndt,nout,nave,num_ndt,print_ndt
       logical decay_p,prod_p,slab_p,cyl_p,sph_p,esc_clear_p
       logical jw_p, jw_dopp_p
       logical pure_dopp_p
       logical cfr_dopp_p,cfr_lrtz_p,cfr_vgt_p
       logical exact_coh_p, exact_incoh_p
       
       logical save_p, read_p
       character*80 oldfile, newfile

       integer n_phi_gauss, n_r_gauss

       common / rescom0 / n_phi_gauss, n_r_gauss
       common / rescom4 / nra,prra,pa,pac
       common / rescom8 /zmaxdp,dzdp,k0,tvac,
     ^              kv,fv,df,av,pcoll,int_ls,nesc
       common / rescom  / zmax,dz,ci,izmax,nout,nave,num_ndt,
     ^              ikmax,ikbmax,
     ^              print_ndt,ndt,decay_p,prod_p,slab_p,cyl_p,sph_p,
     ^              esc_clear_p, jw_p, jw_dopp_p,save_p, read_p,
     ^              pure_dopp_p,cfr_dopp_p,cfr_lrtz_p,cfr_vgt_p,
     ^              exact_coh_p, exact_incoh_p
      common / rescomc / pathname,
     ^              prodfile,filedat,oldfile,newfile
