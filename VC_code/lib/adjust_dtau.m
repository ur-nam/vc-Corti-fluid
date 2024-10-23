function [MP, Uf] = adjust_dtau(MP,fMP,Uf,Uf_a,Uf_b)

    rel_tol = 1e-2;

    err = max(vecnorm(Uf_a(:,fMP.ind) - Uf_b(:,fMP.ind),2,2));

    if err < rel_tol
        Uf = Uf_a;
    end    

    MP.dtau = 0.8*MP.dtau*sqrt(rel_tol/err);
    fprintf(1,"\n * Updated dtau to: %.2e",MP.dtau(601));
end