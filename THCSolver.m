classdef THCSolver
    % Transient Heat Conduction Solver (Completed 1/26/2024)
    
    % This class calculates the HAMR disk temperature. 
    % The steps to run are
        % 1. Define the class
        % 2. Make changes to the parameters as required
        % 3. Make matrices using the function MakeAllMatrices();
        % 4. Run iterations using RunIterations() function 
        % 5. Enjoy the results


    properties(SetAccess = public)
        Vel = 25;   % Disk Linear Velocity
        PowerInput = 8e-4   % Laser power
        NumIterations = 200 % How many iterations should it take to complete the whole run
        NumPerc = 0.5       % Till what point should the laser run till? 0.5 inplies till the midway point
         
        %Final result matrix
        FinalTemps;
        Finalxcoord;

        %Flags
        PlotFlag = 0;       % Flag if the result needs to be plot or not
        MakeMeshFlag = 0;    % Makes and stores new mesh (do only once in a directory).
    end

    properties (Access = private)
        
        % Property Values that affect matrices
        k = 15;             % Thermal Conductivity
        Cv = 3e6;           % Specific Heat capacity per unit volume (= density*heat capacity)
        LaserDia = 30e-9    % Laser diameter
        MulFac = 32;        % Length of simulation box/Laser diameter
        AmTemp = 300;       % Ambient Temperature

        %Boundary Conditions
        PeakPowerIntensity;
        Cooling = 1e2*ones(6,1);
        AreaInput = 5;

        % Mesh Properties
        Nodes;
        Elements;
        ExPts;
        NumNodes;
        NumEl;

        %Iteration Matrics
        Ktb;
        Ktc;
        Qec;
        Cet;
        Qef;

        %Other matrices
        MFluxElementList = 0;
        NodalResults;
    end

    
    
    methods

        function obj = ReadMesh(obj)
            if ~isfile("TFCMesh.mat") || obj.MakeMeshFlag
               obj.MakeMesh();  
            end

            im = load("TFCMesh.mat");
            obj.Nodes = im.NewPList;
            obj.Elements =im.NewTetras;
            obj.NumNodes = size(obj.Nodes,1);
            obj.NumEl = size(obj.Elements,1);

        end

        function obj = MakeAllMatrices(obj)
            obj = obj.ReadMesh();
            obj = obj.SizeNodes();
            obj = obj.MakeKtb();
            obj = obj.getConvSurfaceHeatFlow();
            obj = obj.getCet();
        end

        function obj = change(obj,Str,Val)

            switch Str
                case 'k'
                    obj.k = Val;
                    obj = obj.MakeKtb();

                case 'Cv'
                    obj.Cv = Val;
                    obj = obj.getCet();

                case 'MulFac'
                    obj.MulFac = Val;
                    obj = obj.MakeAllMatrices();

                case 'T0'
                    obj.AmTemp = Val;
                    obj = obj.getConvSurfaceHeatFlow();

                case 'LaserDia'
                    obj.LaserDia = Val;
                    obj = MakeAllMatrices();
                otherwise
                    fprintf('\n Change request unknown');
            end
            
        end

        function MakeMesh(~)

            L = 800;
            W = 320;
            H = 60;
            
            Lsz = [60,25,6,2];
            
            x_lims = [0, 1;
                      0, 1;
                      0, 1;
                      0, 1]*L;
            
            y_lims = [0   , 1;
                      0.25, 0.75;
                      0.45, 0.55;]*W;
            y_lims = [y_lims;158,162];
            
            z_lims = [0   , 1;
                      0.4 , 1;
                      0.75, 1;]*H;
            z_lims = [z_lims;56,60];
            
            AllPoints = zeros(10000,3);
            CurPoints = 0;
            CurChecked = 0;
            
            for rc = 1:3
            
                %Creating block points
                BaseP = cell(3,1);
                BaseP{1} = linspace(x_lims(rc,1),x_lims(rc,2),ceil(diff(x_lims(rc,:))/Lsz(rc))+1);
                BaseP{2} = linspace(y_lims(rc,1),y_lims(rc,2),ceil(diff(y_lims(rc,:))/Lsz(rc))+1);
                BaseP{3} = linspace(z_lims(rc,1),z_lims(rc,2),ceil(diff(z_lims(rc,:))/Lsz(rc))+1);
                
                [Xp,Yp,Zp] = meshgrid(BaseP{1},BaseP{2},BaseP{3});
                AllPoints_N = [Xp(:),Yp(:),Zp(:)];
                NumNew = size(AllPoints_N,1);
                
                % Removing old points that overlap in current block
                if CurPoints > 0
                    DelCol = zeros(CurPoints,1);
                    DelC = 0;
                    for PntCntr = CurChecked+1:CurPoints
                
                        CurP = AllPoints(PntCntr,:);
                        if CurP(1) >= x_lims(rc,1) && CurP(1) <= x_lims(rc,2)
                            if CurP(2) >= y_lims(rc,1) && CurP(2) <= y_lims(rc,2)
                                if CurP(3) >= z_lims(rc,1) && CurP(3) <= z_lims(rc,2)
                                    DelC = DelC+1;
                                    DelCol(DelC) = PntCntr;
                                end
                            end
                        end
                    
                    end
                    DelCol(DelCol==0)=[];
                    AllPoints(DelCol,:) = [];
                    CurChecked = size(AllPoints,1);
                end
                AllPoints(CurPoints+1:CurPoints+NumNew,:) = AllPoints_N;
                CurPoints = size(AllPoints,1);
            
            end
            
            % Refine centerline
            sp = Lsz(rc)/2;
            xcenter = 0:sp:L;
            [Xp,Yp,Zp] = meshgrid(xcenter,W/2,H);
            AllPoints_N = [Xp(:),Yp(:),Zp(:)];
            NumNew = size(AllPoints_N,1);
            AllPoints(CurPoints+1:CurPoints+NumNew,:) = AllPoints_N;
            
            
            DT = delaunayTriangulation(AllPoints);
            
            
            % Convert to 10-node tetrahedrons.
            InitPoints = DT.Points;
            InitTetras = DT.ConnectivityList;
            NumTetras = size(InitTetras,1);
            
            LinCons = [1 2;
                      2 3;
                      3 1;
                      1 4;
                      2 4;
                      3 4];     % Line connections
            NumLines = 6;
            
            NewPList = zeros(15000,4);
            NewTetras = zeros(NumTetras,12);
            NumNewP = 0;
            
            % Get tetra connectivity list
            AdjList = getTetraAdj(InitTetras);
            
            CurTetraPoints = zeros(10,3);
            for TetraCnt = 1:NumTetras
                fprintf('\n Dividing tetrahedra %i',TetraCnt);
                NewTetras(TetraCnt,1) = TetraCnt; % Tetrahedron ID;
                EndNodes = InitTetras(TetraCnt,:);
                Pts = InitPoints(EndNodes,:);
                NewTetras(TetraCnt,12) = 1/6*abs(dot(cross(Pts(2,:)-Pts(1,:),Pts(3,:)-Pts(1,:)),Pts(4,:)-Pts(1,:)));
                
                CurTetraPoints(1:4,:) = Pts;
            
                %New Points
                for j = 1:NumLines
                    p1 = CurTetraPoints(LinCons(j,1),:);
                    p2 = CurTetraPoints(LinCons(j,2),:);
                    midp = (p1+p2)/2;
                    CurTetraPoints(j+4,:) = midp; 
                end
            
                % Insert points into global point lists
                for p = 1:10
                    CurP = CurTetraPoints(p,:);
            
                    AdjacentTetras = AdjList(TetraCnt,:);
                    AdjacentTetras(AdjacentTetras==0) = [];
                    AllPastPoints = NewTetras(AdjacentTetras,2:11);
                    AllPastPoints = AllPastPoints(:);
                    [Pcatch, Pindex] = mancheck(CurP, NewPList(AllPastPoints,2:4),'rows'); %Can be replaced with ismember function.
            
                    if Pcatch
                        NewTetras(TetraCnt,p+1) = AllPastPoints(Pindex);
                    else
                        NumNewP = NumNewP+1;
                        NewPList(NumNewP,:) = [NumNewP,CurP];
                        NewTetras(TetraCnt,p+1) = NumNewP;
                    end
                end
            end
            
            NewPList = NewPList(1:NumNewP,:);
            save("TFCMesh.mat","NewPList","NewTetras");
            
            function AdjList = getTetraAdj(ConnList)
            
                NumTetra = size(ConnList,1);
                AdjList = zeros(NumTetra,50);
                
                for tcn = 2:NumTetra
                    fprintf('\n Dividing tetrahedra %i',tcn);
                    CurTetraList = zeros(1,50);
                    AdjCnt = 0;
                    for pcn = 1:4            
                        [rdet, ~]=find(ConnList(1:tcn-1,:)==ConnList(tcn,pcn));
                        rnum = length(rdet);
                        CurTetraList(AdjCnt+1:AdjCnt+rnum) = rdet; % Add i if search for the future ppoints
                        AdjCnt = AdjCnt+rnum;
                    end
            
                    if AdjCnt > 0 
                        CurTetraList = unique(CurTetraList(1:AdjCnt));
                        CurTetraList(CurTetraList==tcn) = [];
                        AdjCnt = length(CurTetraList);
                        AdjList(tcn,1:AdjCnt) = CurTetraList;
                    end
            
                end
            end
            
            function [flag, ind] = mancheck(myp, allp, ~)
            
                ind = 0;
                NumPts = size(allp,1);
            
                myp = repmat(myp,NumPts,1);
                d = allp == myp;
                d = all(d,2);
                flag = any(d);
                if flag
                    ind = find(d,1,"first");
                end
            
            end
        end

        function obj = SizeNodes(obj)
    
            BlockWidth = obj.LaserDia*obj.MulFac;
            RedFac = BlockWidth/max(abs(obj.Nodes(:,3)));  %Reduction Factor
            
            obj.Nodes(:,2:4) = obj.Nodes(:,2:4)*RedFac;
            obj.Elements(:,12) = obj.Elements(:,12)*(RedFac)^3;
            
            obj.ExPts(1) = min(obj.Nodes(:,2));
            obj.ExPts(2) = max(obj.Nodes(:,2));
            obj.ExPts(3) = min(obj.Nodes(:,3));
            obj.ExPts(4) = max(obj.Nodes(:,3));
            obj.ExPts(5) = min(obj.Nodes(:,4));
            obj.ExPts(6) = max(obj.Nodes(:,4));
           
        end

        function obj = MakeKtb(obj)
    
            obj.NumNodes = size(obj.Nodes, 1);
            
            ValCnt = 0;
            DefSize = obj.NumNodes*10;
            Kmax = zeros(DefSize,1);
            II = zeros(DefSize,1);
            JJ = zeros(DefSize,1);
                
            GradMatrix = obj.getNodeGradient();
                
            obj.NumEl = size(obj.Elements,1);
                    
            for Elc = 1:obj.NumEl
                CursN = obj.Elements(Elc, 2:11);
                CurKlocal = obj.k*getKlocal(obj.Nodes(CursN, 2), obj.Nodes(CursN, 3), obj.Nodes(CursN, 4), obj.Elements(Elc, 12));
                for i = 1:10
                    for j = 1:10
                        ValCnt = ValCnt + 1;
                        II(ValCnt) = obj.Elements(Elc,i+1);
                        JJ(ValCnt) =  obj.Elements(Elc,j+1);
                        Kmax(ValCnt) = CurKlocal(i,j);
                    end
                end
                if DefSize < ValCnt + 600
                    DefSize = DefSize + obj.NumNodes*10;
                    II(DefSize) = 0;
                    JJ(DefSize) = 0;
                    Kmax(DefSize) = 0;
                end
            end
            
            
            obj.Ktb = sparse(II(1:ValCnt), JJ(1:ValCnt), Kmax(1:ValCnt));
            fprintf('\n Created Kbt global Matrix');
           
            
            function Klocal = getKlocal(x,y,z, vol)
                J = [x(2) - x(1), x(3)- x(1), x(4) - x(1);
                     y(2) - y(1), y(3)- y(1), y(4) - y(1);
                     z(2) - z(1), z(3)- z(1), z(4) - z(1)]';
                J2 = inv(J);
                
                Klocal = zeros(10,10);
                vals2 =[ 0.58541020	0.13819660	0.13819660	0.25
                        0.13819660	0.58541020	0.13819660	0.25
                        0.13819660	0.13819660	0.58541020	0.25
                        0.13819660	0.13819660	0.13819660	0.25];
                NumPts2 = size(vals2,1);
               
                for Klocali = 1:10
                    for Klocalj = Klocali:10
                        for PCnt2 = 1:NumPts2
                            delN1 = zeros(3,1);
                            delN2 = zeros(3,1);
                            delN1(1) = GradMatrix(PCnt2, Klocali,1);
                            delN1(2) = GradMatrix(PCnt2, Klocali,2);
                            delN1(3) = GradMatrix(PCnt2, Klocali,3);
                            delN2(1) = GradMatrix(PCnt2, Klocalj,1);
                            delN2(2) = GradMatrix(PCnt2, Klocalj,2);
                            delN2(3) = GradMatrix(PCnt2, Klocalj,3);
        
                            sumVal = dot(J2*delN1,J2*delN2);
                            Klocal(Klocali,Klocalj) = Klocal(Klocali,Klocalj) + vals2(PCnt2,4)*sumVal;
                        end
                        Klocal(Klocalj,Klocali) = Klocal(Klocali,Klocalj);
                    end
                end
                Klocal = Klocal*vol;
            end
            
        end

        function obj = getConvSurfaceHeatFlow(obj)
        
            obj.Ktc = zeros(obj.NumNodes, 1);
            obj.Qec = zeros(obj.NumNodes, 1);
            
            BlockDim = [obj.ExPts(2),obj.ExPts(2),obj.ExPts(4),obj.ExPts(4),obj.ExPts(6),obj.ExPts(6)];
            
            for LoadArea = obj.AreaInput

                ElementList = obj.getElements(BlockDim(LoadArea)*rem(LoadArea-1,2), LoadArea);
                Dir = floor(LoadArea/2+0.5);
                xyDir = [2,3,4];
                xyDir(Dir) = [];
                for Elc = 1:size(ElementList,1)
                    CurNodes = ElementList(Elc, 2:4);
                    x = obj.Nodes(CurNodes,xyDir(1));
                    y = obj.Nodes(CurNodes,xyDir(2));
                    LocalTop =  obj.getElementConvFun(x,y,ElementList(Elc, 5));
                    
                    CurNodes = obj.Elements(ElementList(Elc), 2:11);
                    AddVals = LocalTop*obj.Cooling(LoadArea);
                    obj.Ktc(CurNodes) = obj.Ktc(CurNodes) + AddVals;
                    obj.Qec(CurNodes) = obj.Qec(CurNodes) + AddVals*obj.AmTemp;                    
                end
            end

            obj.Ktc = sparse(1:obj.NumNodes, 1: obj.NumNodes, obj.Ktc);
        
        end

        function obj = getCet(obj)
    
            ValCnt = 0;
            DefSize = obj.NumNodes*10;
            Kmax = zeros(DefSize,1);
            II = zeros(DefSize,1);
            JJ = zeros(DefSize,1);
                
            Clocal = zeros(10,10);
            vals =[ 0.58541020	0.13819660	0.13819660	0.25
                    0.13819660	0.58541020	0.13819660	0.25
                    0.13819660	0.13819660	0.58541020	0.25
                    0.13819660	0.13819660	0.13819660	0.25];
        
            for intI = 1:10
                for intJ = 1:10
                    for PCnt = 1:4
        
                        r = vals(PCnt,1);
                        s = vals(PCnt,2);
                        t = vals(PCnt,3);
        
                        N1 = obj.getDNode(r,s,t,intI,0);
                        N2 = obj.getDNode(r,s,t,intJ,0);
        
                        Clocal(intI, intJ) = Clocal(intI, intJ) + vals(PCnt, 4)*N1*N2;
        
                    end
                end
            end
           
            for Elc = 1:obj.NumEl
                CurKlocal = obj.Cv*obj.Elements(Elc, 12)*Clocal;
                for i = 1:10
                    for j = 1:10
                        ValCnt = ValCnt + 1;
                        II(ValCnt) = obj.Elements(Elc,i+1);
                        JJ(ValCnt) =  obj.Elements(Elc,j+1);
                        Kmax(ValCnt) = CurKlocal(i,j);
                    end
                end
                if DefSize < ValCnt + 600
                    DefSize = DefSize + obj.NumNodes*10;
                    II(DefSize) = 0;
                    JJ(DefSize) = 0;
                    Kmax(DefSize) = 0;
                end
            end
            obj.Cet = sparse(II(1:ValCnt), JJ(1:ValCnt), Kmax(1:ValCnt));
            fprintf('\n Created Cet global Matrix');
           
        end

        function obj = getMassfluxVector2(obj,xPerc)

            LoadArea = 6;
            obj.Qef = zeros(size(obj.Nodes, 1), 1);
        
            yPt = obj.ExPts(4)/2;
            xPt = xPerc/100*obj.ExPts(2);
            fprintf('\n Heat spot center at %4.2f nm',xPt*1e9);
              
            if obj.MFluxElementList == 0
                obj.MFluxElementList = obj.getElements(obj.ExPts(LoadArea), LoadArea);
            end
            Dir = floor(LoadArea/2+0.5);
            xyDir = [2,3,4];
            xyDir(Dir) = [];
        
            sig = obj.LaserDia/sqrt(8*log(2));
            PeakA = obj.PowerInput/(pi*obj.LaserDia^2);
                        
            for Elc = 1:length(obj.MFluxElementList)
                %CurNodes = ElementList(Elc, 2:4);
                CurNodes = obj.Elements(obj.MFluxElementList(Elc), 2:11);
                CurNodes(getDelCols(obj.MFluxElementList(Elc, 5))) = [];
                x = obj.Nodes(CurNodes,xyDir(1));
                y = obj.Nodes(CurNodes,xyDir(2));
                       
                TriArea = obj.getTriArea(x,y);
                qEl = zeros(6,1);
                
                for ri = 1:6
                    radDist = sqrt((x(ri)-xPt)^2 + (y(ri)-yPt)^2);
        
                    % Gaussian Heatflux
                    if radDist < 6*obj.LaserDia
                        qEl(ri) = PeakA*exp(-radDist^2/(2*sig^2));
                    end
                    
                    % Uniform circular heat flux
        %             if radDist < LaserRad
        %                 qEl(ri) = HeatFlux;
        %             end
                end
        
                LocalTop = IntegrateLocalTop(qEl)*TriArea;
                CurNodes = obj.Elements(obj.MFluxElementList(Elc), 2:11);
                CurNodes(getDelCols(obj.MFluxElementList(Elc, 5))) = [];
                obj.Qef(CurNodes) = obj.Qef(CurNodes) + LocalTop;               
                
            end
            
            function DelCol = getDelCols(n)      
            % This function is to calculate the location 
            % of the nodes which are off the surface plane
                switch n
                    case 1
                        DelCol = [1,5,7,8];
                    case 2
                        DelCol = [2,5,6,9];
                    case 3
                        DelCol = [3,6,7,10];
                    case 4
                        DelCol = [4,8,9,10];
                end
            end
        
            function LocalTop = IntegrateLocalTop(qEl)
            % This functin does the integration of the mass heat flux local vector
            
                Vals = [0.8168476 0.0915762 0.0915762 0.1099517;
                        0.0915762 0.8168476 0.0915762 0.1099517;
                        0.0915762 0.0915762 0.8168476 0.1099517;
                        0.1081030 0.4459485 0.4459485 0.2233816;
                        0.4459485 0.1081030 0.4459485 0.2233816;
                        0.4459485 0.4459485 0.1081030 0.2233816];
                    
                Nc = [2 3 4 6 9 10];
                
                qElmid = zeros(7,1);
                for Pindex = 1:6
                   for Hindex = 1:6
                       qElmid(Pindex) = qElmid(Pindex) + qEl(Hindex)*obj.getDNode(Vals(Pindex, 1), Vals(Pindex, 2),Vals(Pindex, 3), Nc(Hindex), 0);
                   end
                end
                    
                LocalTop = zeros(6,1);
                for ci = 1:6
                    for Inti = 1:6
                        LocalTop(ci) = LocalTop(ci) + Vals(Inti, 4)*qElmid(Inti)*obj.getDNode(Vals(Inti, 1), Vals(Inti, 2),Vals(Inti, 3), Nc(ci), 0);
                    end
                end
                
            end
            
        end

        function obj = RunIterations(obj)

            obj.NodalResults =  300*ones(size(obj.Qec));    
            TimeCover = obj.ExPts(2)/obj.Vel;
            dt = TimeCover/obj.NumIterations;
            TotalIters = floor(obj.NumIterations*obj.NumPerc);

            FN = [obj.getNodes(obj.ExPts(1),1); obj.getNodes(obj.ExPts(2),1); obj.getNodes(obj.ExPts(3),2); obj.getNodes(obj.ExPts(4),2)];
            FN = unique(FN);
            uFN = setdiff(obj.Nodes(:,1),FN);

            for i = 1:TotalIters
                obj = obj.getMassfluxVector2(i*obj.Vel*dt/obj.ExPts(2)*100);

                LoadVec = (obj.Qec + obj.Qef)*dt + obj.Cet*obj.NodalResults;
                Kmat = (obj.Ktb + obj.Ktc)*dt + obj.Cet;

                D = Kmat(uFN,uFN);
                C = Kmat(uFN,FN);
                L1 = LoadVec(uFN);
                u = obj.AmTemp*ones(size(FN));
                [v,~,~,~]  = pcg(D,(L1-C*u),1e-6,8000, [],[],obj.NodalResults(uFN));

                obj.NodalResults(uFN) = v;
                obj.NodalResults(FN) = u;

                [~, ~] = obj.PlotTopline();

            end
        end

        function [Xq, Vq] = PlotTopline(obj)
            
            CurrentNodes = obj.getNodes(obj.ExPts(6), 3);
            WantedNodes = obj.Nodes(CurrentNodes, 2:4);
            
            X = WantedNodes(:, 1);
            Y = WantedNodes(:, 2);
            Z = obj.NodalResults(CurrentNodes);
            
            Xq = linspace(0,obj.ExPts(2),500);
            Yq = obj.ExPts(4)/2;
            [Xq, Yq] = meshgrid(Xq,Yq);            
            Vq = griddata(X,Y,Z,Xq,Yq);
            
            if obj.PlotFlag
                plot(Xq,Vq,'LineWidth',2);
                ylim([200 800]);
                pause(0.01);
            end

        end

    end

    methods (Access = private)
        function GradMatrix = getNodeGradient(obj)
            vals =[ 0.58541020	0.13819660	0.13819660	0.25
                    0.13819660	0.58541020	0.13819660	0.25
                    0.13819660	0.13819660	0.58541020	0.25
                    0.13819660	0.13819660	0.13819660	0.25];
        
            NumPts = size(vals,1);
            GradMatrix = zeros(NumPts, 10, 3);
            for i1 = 1:10
                for PCnt = 1:NumPts
                    GradMatrix(PCnt, i1, 1) = obj.getDNode(vals(PCnt,1),vals(PCnt,2),vals(PCnt,3),i1,1);
                    GradMatrix(PCnt, i1, 2) = obj.getDNode(vals(PCnt,1),vals(PCnt,2),vals(PCnt,3),i1,2);
                    GradMatrix(PCnt, i1, 3) = obj.getDNode(vals(PCnt,1),vals(PCnt,2),vals(PCnt,3),i1,3);
                end
            end
        end

        function res = getDNode(~,r,s,t,Ni,Dj)
           
            L1 = 1-r-s-t;
            L2 = r;
            L3 = s;
            L4 = t;
            
            switch Dj
                case 0
                    switch Ni
                        case 1
                            res = (2*L1 - 1)*L1;
                        case 2
                            res = (2*L2 - 1)*L2;
                        case 3
                            res = (2*L3 - 1)*L3;
                        case 4
                            res = (2*L4 - 1)*L4;
                        case 5
                            res = 4*L1*L2;
                        case 6
                            res = 4*L2*L3;
                        case 7
                            res = 4*L1*L3;
                        case 8
                            res = 4*L1*L4;
                        case 9
                            res = 4*L2*L4;
                        case 10
                            res = 4*L3*L4;
                    end
                case 1
                    switch Ni
                        case 1
                            res = 1-4*L1;
                        case 2
                            res = 4*L2-1;
                        case 3
                            res = 0;
                        case 4
                            res = 0;
                        case 5
                            res = 4*(L1-L2);
                        case 6
                            res = 4*L3;
                        case 7
                            res = -4*L3;
                        case 8
                            res = -4*L4;
                        case 9
                            res = 4*L4;
                        case 10
                            res = 0;
                    end                
                case 2
                    switch Ni
                        case 1
                            res = 1-4*L1;
                        case 2
                            res = 0;
                        case 3
                            res = 4*L3-1;
                        case 4
                            res = 0;
                        case 5
                            res = -4*L2;
                        case 6
                            res = 4*L2;
                        case 7
                            res = 4*(L1-L3);
                        case 8
                            res = -4*L4;
                        case 9
                            res = 0;
                        case 10
                            res = 4*L4;
                    end
                case 3
                    switch Ni
                        case 1
                            res = 1-4*L1;
                        case 2
                            res = 0;
                        case 3
                            res = 0;
                        case 4
                            res = 4*L4-1;
                        case 5
                            res = -4*L2;
                        case 6
                            res = 0;
                        case 7
                            res = -4*L3;
                        case 8
                            res = 4*(L1-L4);
                        case 9
                            res = 4*L2;
                        case 10
                            res = 4*L3;
                    end
            end
        end

        function res = getElements(obj, pts, Ar)
        
            Ar = floor((Ar+1)/2);
            NodeList = obj.getNodes(pts, Ar);
            obj.NumEl = size(obj.Elements, 1);
            
            res = zeros(obj.NumEl, 5);
            resC = 0;
            for i = 1:obj.NumEl        
                [c,ia,~] = intersect(obj.Elements(i,2:5),NodeList);
                if length(c) == 3
                    Lones = 1:4;
                    Lones(ia) = [];
                    resC = resC + 1;
                    res(resC, 1) = obj.Elements(i,1);
                    res(resC, 2:4) = c;
                    res(resC, 5) = Lones;
                end
            end
            res = res(1:resC, :);
        
        end

        function res = getNodes(obj, pts, Ar)
    
            obj.NumNodes = length(obj.Nodes);
            res = [];
            for i = 1:obj.NumNodes
                if obj.Nodes(i,Ar+1) == pts
                    res = [res; obj.Nodes(i,1)];
                end
            end
        end

        function Flocal = getElementConvFun(obj,x,y,skip)
            
            J = [x(2) - x(1), x(3)- x(1)
                 y(2) - y(1), y(3)- y(1)]';
            
            if skip == 0
                val = abs(det(J))/12;
                Flocal = val*ones(3,1);
                return;
            end
             
             
            Flocal = zeros(10,1);
            vals = [0.816847573 0.091576214 0.091576214 0.109951744;
                    0.091576214 0.816847573 0.091576214 0.109951744;
                    0.091576214 0.091576214 0.816847573 0.109951744;
                    0.108103018 0.445948491 0.445948491 0.223381590;
                    0.445948491 0.108103018 0.445948491 0.223381590;
                    0.445948491 0.445948491 0.108103018 0.223381590];
            
            switch skip
                case 2
                    vals = [0*vals(:,1), vals(:,2), vals(:,3), vals(:,4)];
                case 3
                    vals = [vals(:,1), 0*vals(:,2), vals(:,3), vals(:,4)];
                case 4
                    vals = [vals(:,1), vals(:,2), 0*vals(:,3), vals(:,4)];
            end
            
            NumPts = size(vals,1);
            
            for Flocali = 1:10
                for PCnt = 1:NumPts
                    sumVal = vals(PCnt,4)*obj.getDNode(vals(PCnt,1),vals(PCnt,2),vals(PCnt,3),Flocali,0);
                    Flocal(Flocali) = Flocal(Flocali) + sumVal;
                end
            end    
            Flocal = Flocal*abs(det(J))/2;       
        end

        function Area = getTriArea(~, x,y)
    
            TriMat = [x(1), y(1), 1;
                      x(2), y(2), 1;
                      x(3), y(3), 1];
            Area = 0.5*abs(det(TriMat));
        
        end
    end
end

