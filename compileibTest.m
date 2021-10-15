% analysis of compileieb data

addpath('./Utils/COMPlib_r1_1');
names = ["AC", "HE", "REA", "DIS", "WEC", "PAS", "TF", "NN", "HF2D", "TMD", "FS", "ROC"];
valid = 0;
eigs_on_pos = 0;
eigs_out = 0;

for i = 1:12
    for j=1:20
        try
            name = strcat(names(i), string(j));
            [A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib(name);
            disp(name);
            
            eigVals = eigs(A);
            valid = valid + size(eigVals,1);
            for k = 1:size(eigVals)
                if real(eigVals(k)) > 0
                    eigs_on_pos = eigs_on_pos + 1;
                end
                if abs(eigVals(k)) > 1
                    eigs_out = eigs_out + 1;
                end
            end
        catch 
            continue;
        end
    end
end
disp(eigs_on_pos);
disp(eigs_out);
disp(valid);
disp(eigs_on_pos / valid);
disp(eigs_out / valid);
