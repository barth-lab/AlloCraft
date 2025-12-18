% subtract corresponding average value from each of the MI types
MIfilterAvg = zeros(size(MI));
MIfilterAvg(MI<randMIlim) = 0;

for i=1:length(MItable1)
    ind1 = MItable1(i,1); ind2 = MItable1(i,2);
    if(MItable1(i,7)==1)
        if(ind_TM(i))
%             MIfilterAvg(ind1,ind2) = (MI(ind1,ind2)-MIavg_bb_bb_TM)/MIstd_bb_bb_TM;
            MIfilterAvg(ind1,ind2) = (MI(ind1,ind2)-MIavg_bb_bb_TM);
        else
%             MIfilterAvg(ind1,ind2) = (MI(ind1,ind2)-MIavg_bb_bb_loop)/MIstd_bb_bb_loop;
            MIfilterAvg(ind1,ind2) = (MI(ind1,ind2)-MIavg_bb_bb_loop);
        end
    elseif(MItable1(i,7)==2)
        if(ind_TM(i))
%             MIfilterAvg(ind1,ind2) = (MI(ind1,ind2)-MIavg_bb_sc_TM)/MIstd_bb_sc_TM;
            MIfilterAvg(ind1,ind2) = (MI(ind1,ind2)-MIavg_bb_sc_TM);
        else
%             MIfilterAvg(ind1,ind2) = (MI(ind1,ind2)-MIavg_bb_sc_loop)/MIstd_bb_sc_loop;
            MIfilterAvg(ind1,ind2) = (MI(ind1,ind2)-MIavg_bb_sc_loop);
        end
    elseif(MItable1(i,7)==3)
        if(ind_TM(i))
%             MIfilterAvg(ind1,ind2) = (MI(ind1,ind2)-MIavg_sc_sc_TM)/MIstd_sc_sc_TM;
            MIfilterAvg(ind1,ind2) = (MI(ind1,ind2)-MIavg_sc_sc_TM);
        else
%             MIfilterAvg(ind1,ind2) = (MI(ind1,ind2)-MIavg_sc_sc_loop)/MIstd_sc_sc_loop;
            MIfilterAvg(ind1,ind2) = (MI(ind1,ind2)-MIavg_sc_sc_loop);
        end
    end
    MIfilterAvg(ind2,ind1) = MIfilterAvg(ind1,ind2);
end
MIfilterAvg(MIfilterAvg<0) = 0;