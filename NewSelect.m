function DB = NewSelect(P,DB,mu)
    %% Preparing Data  
    index = [];
    for i = 1 : size(P.decs,1)
        dist2 = pdist2(real(P.decs(i,:)),real(DB.decs));
        if min(dist2) > 1e-5
            index =[index,i];
        end
    end
    if isempty(index)
        return; 
    elseif length(index) <= mu
       PopNew = P.decs(index,:);
       DB     = [DB,SOLUTION(PopNew)];
       return; 
    end
    
    PopDec = P.decs(index,:);
    PopObj = P.objs(index,:);ObjMSE = P.objmse(index,:);
    PopCon = P.cons(index,:);ConMSE = P.conmse(index,:);
    
    A2Obj  = DB.objs;
    A2Obj  = unique(A2Obj,'rows');
    
    zmin     = min([A2Obj;PopObj],[],1); zmax = max([A2Obj;PopObj],[],1);
    A2Obj    = (A2Obj - zmin)./max(zmax - zmin,10e-10);
    PopObj   = (PopObj - zmin)./max(zmax - zmin,10e-10);
    ObjMSE   = ObjMSE./(max(zmax - zmin,10e-10).^2);
    
    %% Selection
    global beta NI
    Pindex = true(1,size(PopObj,1));
    while length(find(Pindex==0)) < mu 
        % Diversity Contribution
        Pexist   = find(Pindex==1);
        Angle    = acos(1-pdist2(PopObj(Pindex,:),A2Obj,'cosine'));
        Angle    = sort(Angle,2);
        Angle    = Angle(:,1);
        [~,Rank] = sort(Angle,'descend');
        index    = Rank(1:min([ceil(beta*NI),length(find(Pindex==1))]));
        
        % Convergence and Feasibility Contributions
        L       = Pexist(index);
        PopDec_ = PopDec(L,:);
        PopObj_ = PopObj(L,:);ObjMSE_ = ObjMSE(L,:);
        PopCon_ = PopCon(L,:);ConMSE_ = ConMSE(L,:);
        PoF     = Feasible_Probability(PopCon_,ConMSE_);
        Pro = zeros(1,size(PopDec_,1));
        for i = 1:size(PopDec_,1)
            mean   = PopObj_(i,:) - A2Obj;
            sigma  = sqrt(ObjMSE_(i,:));
            Pro(i) = -sum(prod(normcdf((0-mean)./sigma),2))/size(A2Obj,1).*PoF(i);
        end
        [~,Rank] = sort(Pro);
        DB    = [DB,SOLUTION(PopDec(L(Rank(1)),:))];
        A2Obj = (DB.objs - zmin)./max(zmax - zmin,10e-10);
        Pindex(L(Rank(1))) = 0; 
    end
end

function PoF = Feasible_Probability(PopCon,ConMSE)
    [N,M] = size(PopCon);
    PoF    = ones(N,1);
    for i = 1 : N
        for j = 1 : M
               PoF(i) = PoF(i) * normcdf((0-PopCon(i,j))/sqrt(ConMSE(i,j)));
        end
    end
end
