function diagnose_degenerate(folder)
% Test sensitivity within MATLAB itself, without changing the WEiDA formula.
r=load(fullfile(folder,'eigs_diagnosis.mat'));
indices=[10 size(r.matrices,1)]; names={'ordinary_control','circle_boundary'};
variants={'repeat','start_one_ulp','start_small_alternating', ...
    'matrix_diagonal_one_ulp','matrix_offdiagonal_one_ulp','matrix_diagonal_1e12'};
rows=struct([]);
for k=1:numel(indices)
    i=indices(k); A=squeeze(r.matrices(i,:,:)); start=r.operator_inputs{i}(:,2);
    [V,D,flag]=eigs(A,6,'largestabs','StartVector',start);
    baseline=feature(V,D); projector=V(:,1:2)*V(:,1:2)';
    for j=1:numel(variants)
        B=A; z=start;
        switch variants{j}
            case 'start_one_ulp'
                z(1)=z(1)+eps(z(1));
            case 'start_small_alternating'
                z=z+1e-15*(-1).^(1:numel(z))';
            case 'matrix_diagonal_one_ulp'
                B(1,1)=B(1,1)+eps(B(1,1));
            case 'matrix_offdiagonal_one_ulp'
                B(1,2)=B(1,2)+eps(B(1,2)); B(2,1)=B(1,2);
            case 'matrix_diagonal_1e12'
                B(1,1)=B(1,1)+1e-12;
        end
        [W,E,otherflag]=eigs(B,6,'largestabs','StartVector',z);
        actual=feature(W,E);
        row.group=names{k}; row.variant=variants{j};
        row.runtime=version;
        row.matrix_max_change=max(abs(B(:)-A(:)));
        row.start_max_change=max(abs(z-start));
        row.feature_max_change=max(abs(actual-baseline));
        row.feature_change_allowing_global_flip=min(max(abs(actual-baseline)),max(abs(actual+baseline)));
        row.projector_max_change=max(max(abs(W(:,1:2)*W(:,1:2)'-projector)));
        row.eigenvalue_max_change=max(abs(diag(E)-diag(D)));
        row.baseline_top_gap=abs(D(1,1)-D(2,2));
        row.baseline_eigenpair_residual=max(max(abs(A*V-V*D)));
        row.changed_eigenpair_residual=max(max(abs(B*W-W*E)));
        row.converged=(flag==0 && otherflag==0);
        if isempty(rows)
            rows=row;
        else
            rows(end+1)=row;
        end
    end
end
fid=fopen(fullfile(folder,'degenerate_sensitivity.json'),'w');
guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(rows));
disp(jsonencode(rows));
end

function result=feature(V,D)
result=V*diag(D)/sum(diag(D));
if mean(result>0)>.5 || (mean(result>0)==.5 && sum(result(result>0))>-sum(result(result<0)))
    result=-result;
end
end
