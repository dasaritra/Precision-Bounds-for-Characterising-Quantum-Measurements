function [dfiSpec, dfiTr, Lj] = DQFI_single_para(povm,dpovm)

J = size(povm,1);
dim = size(povm{1},1);
Lj={};

for j=1:J
    [vec,val] = eig(povm{j});
    Lj{j} =0;

    for k1=1:dim
        for k2=1:dim
            vk1= vec(:,k1);
            vk2= vec(:,k2);
            Lj{j} =Lj{j} +2* (vk1*vk1'*dpovm{j} * vk2*vk2')/(val(k1,k1)+ ...
                                                           val(k2,k2));
        end
    end
end
    

dfiMat=0;
dfiTr=0;
for j=1:J
    dfiMat = dfiMat + Lj{j}*povm{j}*Lj{j} ;
    dfiTr= dfiTr + trace(  Lj{j}*povm{j}*Lj{j} );
end


dfiSpec = norm(dfiMat);


    

