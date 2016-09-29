% Elias Cohenca - 302025358
%
% Calculates the Sturm sequence of polynom p between its bounds
% Input: p - polynom (coefficients from high to low degree)
%        bounds - low and high bound to search for roots (optional)
% Output: ret - Roots intervals
function [ret] = SturmRoot(p,bounds)
    tol=1e-10;
    Eps=1e-5; % Interval size
    MaxRecursion=50;
    if nargin == 0
        error('SturmRoot requires arguments');
    end
    p=NormalizePoly(p);
    % Bounds on roots of polynomials
    if exist('bounds','var')
        LowB=bounds(1); HighB=bounds(2);
    else
        LowB=abs(p(end))/(abs(p(end))+max(abs(p(1:end-1))));
        HighB=1+max(abs(p(2:end)))/abs(p(1));
        % This lowerbound is only for roots>0, since these rules are for abs(root)
        if FindSignChangesNeg(p) > 0 || p(end)==0   % Has at least one negative root or zero
            LowB=-HighB;  
        end
    end

    % Main functions
    % Find the sturm sequence for the bounds
    sturm(p,LowB,HighB,false);
    
    retCell=cell(0,1); % init cell for intervals
    % Find the intervals where the roots are starting with the known bounds
    % of the roots
    retCell=FindRootsInterval(p,LowB,HighB,retCell);
    ret=zeros(numel(retCell),1);
    for j=1:numel(retCell)
        ret(j)=(retCell{j}(1)+retCell{j}(2))/2; % Return as double, middle of interval
%         mid=(retCell{j}(1)+retCell{j}(2))/2;
%         ret(j)=NewtonRaphson(p,retCell{j}(1),retCell{j}(2),mid,1);
    end

    % Find the intervals where the roots are.
    % This is a recursive function that divides each interval in two parts,
    % and keeps subdividing the sections where roots are found, until reaching
    % a size of Eps. We know then that in this interval we'll find the root
    function [ret] = FindRootsInterval(p,LowB,HighB,ret)
        delta=HighB-LowB;
        if delta <= Eps
            ret{end+1}=[LowB, HighB];
            return;
        end
        div=2;
        while abs(EvalPoly(p,LowB+delta/div))<tol %root at division limit?
            div=div+1; % Try another division point (1/3, 1/4, etc...)
            if div>100
                ret{end+1}=[LowB, HighB];
                return; % Avoid too small division. Say this is the interval.
            end
        end
        split=LowB+delta/div;
        RootsLeft=sturm(p,LowB,split,false);
        RootsRight=sturm(p,split,HighB,false);
        if RootsLeft > 0
            ret=FindRootsInterval(p,LowB,split,ret);
        end
        if RootsRight > 0
            ret=FindRootsInterval(p,split,HighB,ret);
        end
    end

    % Find the sturm sequence. The delta number of sign changes give the number
    % of roots in the interval.
    % Input: p - the polynom to work on
    %        LowB - the low bound for the sturm sequence
    %        HighB - the high bound for the sturm sequence
    %        printSeq - if true prints the sturm sequence for this interval
    function [ret] = sturm(p,LowB,HighB,printSeq)
        n=numel(p);
        St=zeros(2,n); % To keep the sturm sequence
        seq=zeros(n);
        seq(1,:)=p;
        seq(2,2:n)=PolynomDeriv(p);
        for i=3:n
            first=seq(i-2,find(seq(i-2,:),1,'first'):n); %remove leading zeros
            second=seq(i-1,find(seq(i-1,:),1,'first'):n); %remove leading zeros
            [~,r]=deconv(first,second);
            seq(i,i:n)=-r(end-n+i:end);
        end
        for i=1:n
            St(1,i)=EvalPoly(seq(i,i:end),LowB);
            St(2,i)=EvalPoly(seq(i,i:end),HighB);
        end
        S_Low=FindSignChanges(St(1,:));
        S_High=FindSignChanges(St(2,:));

        ret=S_Low-S_High;
        if printSeq
            str=sprintf('Sturm sequence in interval [%.4f, %.4f]', LowB, HighB);
            disp(str);
            len=length(str);
            disp(repmat('-',1,len));
            for i=1:n
                str=sprintf('P%d\t%12.4f\t|\t%12.4f',i-1,St(1,i),St(2,i));
                disp(str);
            end
            disp(repmat('-',1,len));
            str=sprintf('Sign:%12d\t|\t%12d',S_Low,S_High);
            disp(str);
            str=sprintf('\nNumber of distinct real roots');
            disp(str);
            disp(ret);
        end
    end


    function [newX]=NewtonRaphson(P,a,b,x,rec)
        Q=PolynomDeriv(P); % derivative
        f=EvalPoly(P,x);   % value of function at this point
        der_f=EvalPoly(Q,x);  % value of derivative at this point
        newX=x-f/der_f;
        % Prevent finding a root outside this interval
        if newX>b || newX<a
            newX=x; % Since in this function we deal with very small interval
                    % let's just return the previous best find
            return;
        end

        % If the change is more than 1%, keep going
        if abs((newX-x)/newX)>0.01 && rec<MaxRecursion
            newX=NewtonRaphson(P,a,b,newX,rec+1);
        end
    end

    % Input: p - polynom
    %        tol - tolerance to check for zero
    % Output: sc - sign changes
    function [sc] = FindSignChanges(p)
        sc=numel(find(diff(sign(p(abs(p)>tol)))~=0)); % num of positive roots (or less by 2n)
    end

    % Calculates the first derivative of the polynom p
    function [der] = PolynomDeriv(p)
        n=numel(p);
        der=zeros(1,n-1);
        for i=1:n-1
            der(i)=p(i)*(n-i);
        end
    end

    % Find if has negative roots
    function [sc] = FindSignChangesNeg(p)
        p(end-1:-2:1)=p(end-1:-2:1)*-1;  % Polynom of p(-x)
        sc=FindSignChanges(p);
    end

    % Find value of polynom at x using Horner's rule (more precise than
    % just calculating each coefficient separately)
    function [val] = EvalPoly(p,x)
        val(1:numel(x))=p(1);
        for i=2:numel(p)
            val=val.*x+p(i);
        end
    end

    % Normalize the polynom if needed (the first coefficient != 1)
    function [ret] = NormalizePoly(p)
        if numel(p) > 1 && p(1)~=1
            ret=p/p(1); % Normalize polynom
        else
            ret=p;
        end
    end
end
