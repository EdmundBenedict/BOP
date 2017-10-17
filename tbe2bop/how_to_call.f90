        call sylmnc(cy,16)
        call scg(9,cg,indxcg,jcg)
        call makcg9(indxcg,jcg,cg,gaunt,ak)

        if (.not. allocated(tbc%esamap)) then
            allocate(tbc%esamap(0:nproc), pidx(nbas))
            xl = 1
            if (nlmq==nlmq1) xl = 2
            pidx(1) = (lmxl(s_ctrl%ipc(1))+xl)**2
            do i = 2, nbas
            pidx(i) = pidx(i-1) + (lmxl(s_ctrl%ipc(i))+xl)**2
            end do
            call vbdist(nbas, pidx, nproc, tbc%esamap)
            deallocate(pidx)
            allocate(tbc%escount(0:nproc-1))
            tbc%escount = tbc%esamap(1:nproc) - tbc%esamap(0:nproc-1)
        end if

        allocate(struxidx(nbas, tbc%escount(pid)))
        if (pv) then
            allocate(dstrxidx(nbas,tbc%escount(pid)))
        else
            allocate(dstrxidx(1,1))
        end if
        call mkstrxidx(tbc, pv, force, nbas, s_ctrl%ipc, lmxl, struxsize, dstrxsize, struxidx, dstrxidx)
        memstr = nbas*tbc%escount(pid)
        if (pv) memstr = memstr*2
        memstr = (memstr*4 + (struxsize+dstrxsize)*8)/2**20
        if (memstr > 24) call info2(10,1,1,' TBZINT: allocating  %i MiB memory for strux&dstrx + indices at root',memstr,0)
        allocate(struxd(struxsize))
        allocate(dstrxd(dstrxsize))
        call mkstrxd(s_ctrl,s_ctrl%ipc,s_lat,tbc,nlmq, nlmq1,lmxl,ldip,dlv2,nkd,qlv2,nkq, &
                     indxcg,jcg,cg,cy,struxd,dstrxd, struxidx, dstrxidx, struxsize, dstrxsize)



