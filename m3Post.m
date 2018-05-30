%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% PostProcess  File.
%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function  [resultList] = M3_addZeroZ(input_list)
%
% FUNCTION NOTE: 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

[num samples] = size(input_list);

%resultList = zeros(1,samples);
resultList = sum(input_list, 1);

for i=1:num
    idx  = input_list(i,:) < e^-9;
    resultList(idx) = 0;
end

return;
endfunction






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function  [resultList] = M3_addZero(input_list1, input_list2)
%
% FUNCTION NOTE: 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

len1= length(input_list1);
len2= length(input_list2);
resultList = zeros(1,len1);

if (len1 != len2)
    return;
end

resultList = input_list1 + input_list2;

idx1 = input_list1< e^-9;
idx2 = input_list2< e^-9;

resultList(idx1) = 0;
resultList(idx2) = 0;
return;
endfunction




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function  [resultList] = M3_regular(input_list)
%
% FUNCTION NOTE: 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

resultList = input_list;
tMax = max(input_list);
tMin = min(input_list);
gg = tMax - tMin;
if (gg<=0)
    return;
end
%resultList = input_list ./ tMax;
resultList = (input_list - tMin) ./ gg;

return;
endfunction





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function  [resultList] = M3_rectangle(input_list, greater, numbers,     negative)
%
% FUNCTION NOTE: 
% greater > 0  means select the greater ones.
%
% Output:
%  0 --- not selected.
%  1 --- selected.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

len = length(input_list);
selectIdx=0;

if (greater > 0)
    [a selectIdx] = sort(input_list, "descend");
else
    [a selectIdx] = sort(input_list);
end

if (negative > 0)
    resultList = zeros(size(input_list));
    resultIdx = selectIdx(1:numbers);
    resultList(resultIdx) = 1;
else
    resultList = ones(size(input_list));
    %resultIdx = selectIdx((numbers+1):len);
    resultIdx = selectIdx(1:numbers);
    resultList(resultIdx) = 0;
end


return;
endfunction



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function  [resultList] = M3_rectangle_mannul(input_list)
%
% FUNCTION NOTE: 
% greater > 0  means select the greater ones.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

len = length(input_list);

resultList = zeros(1,len);
valididx = input_list>0;
resultList(valididx) = 1;

return;
endfunction






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function  [resultIdx] = M3_getFilter(input_list, greater, numbers,     negative)
%
% FUNCTION NOTE: 
% greater > 0  means select the greater ones.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

len = length(input_list);
selectIdx=0;

if (greater > 0)
    [a selectIdx] = sort(input_list, "descend");
else
    [a selectIdx] = sort(input_list);
end

if (negative > 0)
    resultIdx = selectIdx((numbers+1):len);
else
    resultIdx = selectIdx(1:numbers);
end

return;
endfunction



######################################################## FILE.
function  [hcl] = M3_getHcl(input_list)
######################################################## FILE.


hcl = 0;
list = input_list;
rawListLen = length(list);
pastMax(1) = list(1);

if (sum(list<e^-6) > 0)
    return;
end


for i=2:rawListLen
    v = list(i);
    maxThis = pastMax(i-1);
    if (v>maxThis)
        maxThis = v;
    end
    pastMax(i) = maxThis;
end

%input_list
%curHcl = (list'-pastMax) ./ pastMax
%hcl = 1 * min(curHcl)
curHcl = (list'-pastMax) ./ pastMax;
hcl = min(curHcl);
return;

endfunction


######################################################## FILE.
function  [liOut] = M3_getAver(input_list, interval)
######################################################## FILE.


liOut = 0;
list = input_list;
rawListLen = length(list);

for i=1:rawListLen
    startIdx = i-interval;
    endIdx   = i;
    startIdx = max(1, startIdx);
    liOut(i) = mean(list(startIdx:endIdx));
end

return;
endfunction


######################################################## FILE.
function  [activityCount] = M3_getActivity(list, avList)
######################################################## FILE.


activityCount = 0;
rawListLen = length(list);
rawListLen2= length(avList);
if (rawListLen != rawListLen2)
    return ;
end


count = 0;
state = 0;
for i=1:rawListLen
    if (state > 0) 
        if (list(i) < avList(i) )
            state = -1;
            count++;
        end
    else
        if (list(i) > avList(i) )
            state = 1;
            count++;
        end
    end 
end

activityCount = count;
return;

endfunction




######################################################## FILE.
function  [output_roi  output_time] = M3_getRoiOfFuture_avline(input_ap, input_sid)
######################################################## FILE.
global m3 ;
global m3_ids ;
global nGroups;

global m3_header_width;
global m3_index_date1     ;    % header area
global m3_index_date2     ;    % header area
global m3_index_market    ;   % header area
global m3_index_strategy2_y_estimate;
global m3_index_strategy2_y_estimate2;
global m3_id_width;
global m3_index_strategy1_roi;

global m3_index_sid ; % sid 
global stockIdStr;

global m2;
global m2_header_width;
global m2_index_p;
global m2_index_v;
global m2_id_width;


thisRoi = 0;
output_roi = 0;
output_time = 0;
sid = input_sid;
end_ap = input_ap + 5;

indexRoi = m3_header_width + m3_index_strategy1_roi + (sid-1)*m3_id_width;
m3(input_ap, indexRoi) = 0; % initialize to 0.

if (end_ap  > nGroups)
    end_ap = nGroups;
    %return;
end


rawListLen = 0;
rawList = 0;
for ap=input_ap:end_ap
       % ******************************************************
       %   generate roi, and store to m3_index_strategy1_roi
       % ******************************************************
       thisAPidx  = (m2(:,2) == ap);
       rawPIdx    =  m2_header_width + m2_index_p - 2 + (sid-1)*m2_id_width;
       rawP       =  m2(thisAPidx, rawPIdx);
       samples = length(rawP);
       quaterN = samples/4;

       if (rawListLen < 1)
           rawList =           rawP(1 : quaterN*2);
       else
           rawList = [rawList; rawP(1 : quaterN*2)];
       end
       rawListLen = length(rawList);
end


if (rawListLen<100)
    rawListLen
    return;
end

avList     = M3_getAver(rawList, 150);
rawList    = rawList(101:rawListLen);
avList     = avList(101:rawListLen);
rawListLen = rawListLen - 100;


input_ap
sid
rawListLen
if (rawListLen<10)
    return;
end
vv_first= rawList(1)
if (vv_first < e^-10)
    return;
end


vv_zhisun = vv_first * (1-0.03)
vv_zhiyin = vv_first * (1+0.03)
vv_buyPos = vv_first * (1-0.003)
lastV  = vv_first;
for i=28:(rawListLen-0)
       v = rawList(i);

       output_time = i;

       if (v < e^-5)
           thisRoi = 0;
           break;
       end
       if (abs((v - lastV)/v) > 0.18)
           thisRoi = 0;
           break;
       end
       lastV = v;


       if (i<30)
           continue;
       end

       
       thisRoi = (v-vv_first)/vv_first;
       %hcl = M3_getHcl(rawList(30:i))

       if (v < avList(i))
           breakCause = 77
           breakHere = i
           break;
       end

       % add at 2018.01.13
       %if (i==51 && v>vv_first)
           %breakCause = 3
           %breakHere = i
           %break;
       %end

       %if (v>vv_zhiyin)
           %breakCause = 4
           %breakHere = i
           %break;
       %end

       %if (v<vv_zhisun)
           %breakCause = 5
           %breakHere = i
           %break;
       %end

       %if (thisRoi>0 && hcl < -0.02)
       %if (hcl < -0.02)
           %breakCause = 5
           %breakHere = i
           %break;
       %end

       %if (thisRoi<0 && hcl < -0.01)
           %breakCause = 7
           %breakHere = i
           %break;
       %end


       %if (i==82)
           %breakCause = 5
           %breakHere = i
           %break;
       %end
end


thisRoi


% ******************************************************
%   store to m3_index_strategy1_roi
% ******************************************************
indexRoi = m3_header_width + m3_index_strategy1_roi + (sid-1)*m3_id_width;
m3(input_ap, indexRoi) = thisRoi;


% ******************************************************
%   return 
% ******************************************************
output_roi = thisRoi;
output_time ;
return;

endfunction





######################################################## FILE.
function  [output_roi  output_time] = M3_getRoiOfFuture_dynSell_avg(input_ap, input_sid)
######################################################## FILE.
global m3 ;
global m3_ids ;
global nGroups;

global m3_header_width;
global m3_index_date1     ;    % header area
global m3_index_date2     ;    % header area
global m3_index_market    ;   % header area
global m3_index_strategy2_y_estimate;
global m3_index_strategy2_y_estimate2;
global m3_id_width;
global m3_index_strategy1_roi;

global m3_index_sid ; % sid 
global stockIdStr;

global m2;
global m2_header_width;
global m2_index_p;
global m2_index_v;
global m2_id_width;


thisRoi = 0;
output_roi = 0;
output_time = 0;
sid = input_sid;
end_ap = input_ap + 1;
%end_ap = input_ap + 5; % lazy policy: very dangerous.

indexRoi = m3_header_width + m3_index_strategy1_roi + (sid-1)*m3_id_width;
m3(input_ap, indexRoi) = 0; % initialize to 0.

if (end_ap  > nGroups)
    end_ap = nGroups;
    %return;
end


rawListLen = 0;
rawList = 0;
for ap=input_ap+1:end_ap

       % ******************************************************
       %   generate roi, and store to m3_index_strategy1_roi
       % ******************************************************
       thisAPidx  = (m2(:,2) == ap);
       rawPIdx    =  m2_header_width + m2_index_p - 2 + (sid-1)*m2_id_width;
       rawP       =  m2(thisAPidx, rawPIdx);
       samples = length(rawP);
       quaterN = samples/4;

       if (rawListLen < 1)
           rawList =           rawP(quaterN*2+1 : quaterN*4);
       else
           rawList = [rawList; rawP(quaterN*2+1 : quaterN*4)];
       end
       rawListLen = length(rawList);
end



input_ap
sid
rawList
rawListLen
if (rawListLen<10)
    return;
end
vv_first= rawList(1)
%vv_first= rawList(10)
if (vv_first < e^-10)
    return;
end


vv_zhisun = vv_first * (1-0.03)
vv_zhiyin = vv_first * (1+0.03)
vv_buyPos = vv_first * (1-0.003)
lastV  = vv_first;
successBuy = 0;
successBuy = 1;
keyRoi = 0;
if (1)
for i=30:(rawListLen-0)
       v = rawList(i);
       output_time = i;

       if (v < e^-5)
           thisRoi = 0;
           break;
       end
       if (abs((v - lastV)/v) > 0.18)
           thisRoi = 0;
           break;
       end

       lastV = v;
       thisRoi = (v-vv_first)/vv_first;
if (1)
       hcl = M3_getHcl(rawList(30:i))

       if (i == 30)
           keyRoi = thisRoi;
       end
       
       if (keyRoi < 0.01 && thisRoi > keyRoi+0.01)
           %breakCause = 6
           %breakHere = i
           %break;
       end

       if (hcl<-0.02)
           breakCause = 7
           breakHere = i
           break;
       end
end


       % add at 2018.01.13
       %if (i==51 && v>vv_first)
           %breakCause = 3
           %breakHere = i
           %break;
       %end
       %if (v>vv_zhiyin)
           %breakCause = 4
           %breakHere = i
           %break;
       %end
       %if (v<vv_zhisun)
           %breakCause = 5
           %breakHere = i
           %break;
       %end
       %if (thisRoi>0 && hcl < -0.02)
       %if (hcl < -0.02)
           %breakCause = 5
           %breakHere = i
           %break;
       %end
       %if (thisRoi<0 && hcl < -0.01)
           %breakCause = 7
           %breakHere = i
           %break;
       %end
       %if (i==82)
           %breakCause = 5
           %breakHere = i
           %break;
       %end
end
end


thisRoi


% ******************************************************
%   store to m3_index_strategy1_roi
% ******************************************************
indexRoi = m3_header_width + m3_index_strategy1_roi + (sid-1)*m3_id_width;
m3(input_ap, indexRoi) = thisRoi;


% ******************************************************
%   return 
% ******************************************************
output_roi = thisRoi;
output_time ;
return;

endfunction




######################################################## FILE.
function  [output_roi  output_time] = M3_getRoiOfFuture_dynSell_avg_afternoon(input_ap, input_sid)
######################################################## FILE.
global m3 ;
global m3_ids ;
global nGroups;

global m3_header_width;
global m3_index_date1     ;    % header area
global m3_index_date2     ;    % header area
global m3_index_market    ;   % header area
global m3_index_strategy2_y_estimate;
global m3_index_strategy2_y_estimate2;
global m3_id_width;
global m3_index_strategy1_roi;

global m3_index_sid ; % sid 
global stockIdStr;

global m2;
global m2_header_width;
global m2_index_p;
global m2_index_v;
global m2_id_width;


thisRoi = 0;
output_roi = 0;
output_time = 0;
sid = input_sid;
end_ap = input_ap + 1;
%end_ap = input_ap + 5; % lazy policy: very dangerous.

indexRoi = m3_header_width + m3_index_strategy1_roi + (sid-1)*m3_id_width;
m3(input_ap, indexRoi) = 0; % initialize to 0.

if (end_ap  > nGroups)
    end_ap = nGroups;
    %return;
end


rawListLen = 0;
rawList = 0;
ap=input_ap

       % ******************************************************
       %   generate roi, and store to m3_index_strategy1_roi
       % ******************************************************
       thisAPidx  = (m2(:,2) == ap);
       rawPIdx    =  m2_header_width + m2_index_p - 2 + (sid-1)*m2_id_width;
       rawP       =  m2(thisAPidx, rawPIdx);
       samples = length(rawP);
       quaterN = samples/4;

       if (rawListLen < 1)
           rawList =           rawP(quaterN*2+1 : quaterN*4);
       else
           rawList = [rawList; rawP(quaterN*2+1 : quaterN*4)];
       end
       rawListLen = length(rawList);
%end



input_ap
sid
rawList
rawListLen
if (rawListLen<10)
    return;
end
%vv_first= rawList(24)
%vv_first= rawList(10)
vv_first= rawList(13)
if (vv_first < e^-10)
    return;
end


vv_zhisun = vv_first * (1-0.03)
vv_zhiyin = vv_first * (1+0.03)
lastV  = vv_first;
if (1)
%for i=1:rawListLen
for i=30:(rawListLen-0)
       v = rawList(i);

       output_time = i;

       if (v < e^-5)
           thisRoi = 0;
           break;
       end
       if (abs((v - lastV)/v) > 0.18)
           thisRoi = 0;
           break;
       end
       lastV = v;

       
       thisRoi = (v-vv_first)/vv_first;

       if (1)
       dur_start = i - 20;
       dur_start = i - 5;
       dur_start = i - 2;
       dur_end   = i ;
       %dur_end   = i - 5;
       %vVS = 0.98*mean(rawList(dur_start:dur_end)) ;
       %vVS = 1.025*mean(rawList(dur_start:dur_end)) ;
       %vVS = 1.015*mean(rawList(dur_start:dur_end)) ;
       %vVS = 1.025*mean(rawList(dur_start:dur_end)) ;
       %vVS = 1.05*mean(rawList(dur_start:dur_end)) ;
       vVS = rawList(dur_start) * 1.025;
       if (v > vVS)
           %breakCause = 2
           %breakHere = i
           %break;
       end

       if (i==40)
           %breakCause = 7
           %breakHere = i
           %break;
       end



       % add at 2018.01.13
       if (i==51 && v>vv_first)
           %breakCause = 3
           %breakHere = i
           %break;
       end

       if (v>vv_zhiyin)
           %breakCause = 4
           %breakHere = i
           %break;
       end

       if (v<vv_zhisun)
           %breakCause = 5
           %breakHere = i
           %break;
       end

       if (i==82)
           %breakCause = 5
           %breakHere = i
           %break;
       end
       end
end
end


thisRoi


% ******************************************************
%   store to m3_index_strategy1_roi
% ******************************************************
indexRoi = m3_header_width + m3_index_strategy1_roi + (sid-1)*m3_id_width;
m3(input_ap, indexRoi) = thisRoi;


% ******************************************************
%   return 
% ******************************************************
output_roi = thisRoi;
output_time ;
return;

endfunction







######################################################## FILE.
function  [output_roi] = M3_getRoiOfFuture_dynSell(input_ap, input_sid)
######################################################## FILE.
global m3 ;
global m3_ids ;
global nGroups;

global m3_header_width;
global m3_index_date1     ;    % header area
global m3_index_date2     ;    % header area
global m3_index_market    ;   % header area
global m3_index_strategy2_y_estimate;
global m3_index_strategy2_y_estimate2;
global m3_id_width;
global m3_index_strategy1_roi;

global m3_index_sid ; % sid 
global stockIdStr;

global m2;
global m2_header_width;
global m2_index_p;
global m2_index_v;
global m2_id_width;


thisRoi = 0;
output_roi = 0;
sid = input_sid;
%end_ap = input_ap + 1;
end_ap = input_ap + 5; % lazy policy: very dangerous.

indexRoi = m3_header_width + m3_index_strategy1_roi + (sid-1)*m3_id_width;
m3(input_ap, indexRoi) = 0; % initialize to 0.

if (end_ap  > nGroups)
    end_ap = nGroups;
    %return;
end


rawListLen = 0;
rawList = 0;
for ap=input_ap+1:end_ap

       % ******************************************************
       %   generate roi, and store to m3_index_strategy1_roi
       % ******************************************************
       thisAPidx  = (m2(:,2) == ap);
       rawPIdx    =  m2_header_width + m2_index_p - 2 + (sid-1)*m2_id_width;
       rawP       =  m2(thisAPidx, rawPIdx);
       samples = length(rawP);
       quaterN = samples/4;

       if (rawListLen < 1)
           rawList =           rawP(quaterN*2+1 : quaterN*4);
       else
           rawList = [rawList; rawP(quaterN*2+1 : quaterN*4)];
       end
       rawListLen = length(rawList);
end



input_ap
sid
rawList
rawListLen
vv_first= rawList(1)
if (vv_first < e^-10)
    return;
end
maxHcl  = -0.02;


lastV  = vv_first;
tmpMax = 0;
for i=1:rawListLen
       v = rawList(i);

       if (v < e^-5)
           thisRoi = 0;
           break;
       end
       if (abs((v - lastV)/v) > 0.18)
           thisRoi = 0;
           break;
       end
       lastV = v;

       %vZL = v - vv_first*(1+(0.04/50)*i);
       vZL = v - vv_first*(1+(0.02/50)*i);
       if (vZL>tmpMax)
           tmpMax = vZL;
       end
       %hcl = (vZL-tmpMax) / tmpMax;
       hcl = (vZL-tmpMax) / vv_first;
       
       thisRoi = (v-vv_first)/vv_first;

       if (hcl < maxHcl)
           hcl
           breakHere = i
           break;
       end

end


thisRoi


% ******************************************************
%   store to m3_index_strategy1_roi
% ******************************************************
indexRoi = m3_header_width + m3_index_strategy1_roi + (sid-1)*m3_id_width;
m3(input_ap, indexRoi) = thisRoi;


% ******************************************************
%   return 
% ******************************************************
output_roi = thisRoi;
return;

endfunction









######################################################## FILE.
function  [output_roi] = M3_getRoiOfFuture(input_ap, input_sid)
######################################################## FILE.
global m3 ;
global m3_ids ;
global nGroups;

global m3_header_width;
global m3_index_date1     ;    % header area
global m3_index_date2     ;    % header area
global m3_index_market    ;   % header area
global m3_index_strategy2_y_estimate;
global m3_index_strategy2_y_estimate2;
global m3_id_width;
global m3_index_strategy1_roi;

global m3_index_sid ; % sid 
global stockIdStr;

global m2;
global m2_header_width;
global m2_index_p;
global m2_index_v;
global m2_id_width;


thisRoi = 0;
output_roi = 0;
sid = input_sid;
end_ap = input_ap + 1;
%end_ap = input_ap + 5; % lazy policy: very dangerous.
%end_ap = input_ap + 2; % also dangerous!

indexRoi = m3_header_width + m3_index_strategy1_roi + (sid-1)*m3_id_width;
m3(input_ap, indexRoi) = 0; % initialize to 0.

if (end_ap  > nGroups)
    end_ap = nGroups;
    %return;
end


rawListLen = 0;
rawList = 0;
rawList_key=-1;

for ap=input_ap+1:end_ap

       % ******************************************************
       %   generate roi, and store to m3_index_strategy1_roi
       % ******************************************************
       thisAPidx  = (m2(:,2) == ap);
       rawPIdx    =  m2_header_width + m2_index_p - 2 + (sid-1)*m2_id_width;
       rawP       =  m2(thisAPidx, rawPIdx);
       samples = length(rawP);
       quaterN = samples/4;

       if (rawList_key < 0) 
           rawList_key = [rawP(quaterN*2+1)  rawP(samples)];
       else
           tt_ = [rawList_key  rawP(samples)];
           rawList_key = tt_;
       end


       if (rawListLen < 1)
           rawList = rawP(quaterN*2+1 : quaterN*4);
       else
           rawList = [rawList; rawP(quaterN*2+1 : quaterN*4)];
           %rawList = [rawList; rawP];
       end
       rawListLen = length(rawList);

end



% ******************************************************
%   use the rawList_key to substitude to rawList.
%    IIR mechanism.
% ******************************************************
if (1)
    rawList = rawList_key 
    rawListLen = length(rawList)
end




input_ap
sid
rawList
vv_first= rawList(1)
if (vv_first < e^-10)
    return;
end
vv_up   = vv_first * (1 + 0.08);
vv_down = vv_first * (1 - 0.05);


lastV = vv_first;
for i=1:rawListLen

       v = rawList(i);

       if (v < e^-5)
           thisRoi = 0;
           break;
       end
       if (abs((v - lastV)/v) > 0.18)
           thisRoi = 0;
           break;
       end
       lastV = v;
       

       thisRoi = (v-vv_first)/vv_first;

       if (v>= vv_up)
           %break;
       end
       if (v<= vv_down)
           %break;
       end

       if (i>1)
       if (rawList(i) < rawList(i-1))
           breakHere = i
           break;
       end
       end

       % ******************************************************
       %   if no break, then fixed sell point.
       % ******************************************************
end


thisRoi

% ******************************************************
%   store to m3_index_strategy1_roi
% ******************************************************
indexRoi = m3_header_width + m3_index_strategy1_roi + (sid-1)*m3_id_width;
m3(input_ap, indexRoi) = thisRoi;


% ******************************************************
%   return 
% ******************************************************
output_roi = thisRoi;
return;

endfunction




######################################################## FILE.
function  [] = M3_compositeSids()
######################################################## FILE.
global m3 ;
global m3_ids ;
global nGroups;

global m3_header_width;
global m3_index_date1     ;    % header area
global m3_index_date2     ;    % header area
global m3_index_market    ;   % header area
global m3_index_strategy2_y_estimate;
global m3_index_strategy2_y_estimate2;
global m3_id_width;

global m3_index_sid ; % sid 
global stockIdStr;
global stockBriefIDStr;


for ap=nGroups-1:nGroups
    indexStart =  m3_header_width + m3_index_sid;
    indexEnd   =  m3_header_width + m3_ids*m3_id_width;
    Idx       =  indexStart : m3_id_width : indexEnd; 

    % *******************
    %  fill.
    % *******************
    stockIdStr = m3(ap, Idx);
end


for sid=1:m3_ids
    t = stockBriefIDStr(sid);
    if (t<300000)
        stockBriefIDStr(sid) = t + 900000;
    end
end

if (stockBriefIDStr != stockIdStr')
    exceptionErrUnMatch = 0
    exit 0
end


endfunction

######################################################## FILE.
function  [] = M3_postRank()
######################################################## FILE.
global m2 ;
global m3 ;
global m3_ids ;
global nGroups;

global m3_header_width;
global m3_index_date1     ;    % header area
global m3_index_date2     ;    % header area
global m3_index_market    ;   % header area
global m3_index_market_1  ;
global m3_index_market_2  ; 
global m3_index_market_3  ; 
global m3_index_strategy2_y_estimate;
global m3_index_strategy2_y_estimate2;
global m3_id_width;

global m3_index_strategy2_b0;
global m3_index_strategy2_b0_pvalue ;
global m3_index_strategy2_y ;
global m3_index_strategy1_roi;
global m3_index_strategy1_y ;
global m3_index_strategy2_indicator1 ;

global m3_index_global_valid;

global m2_header_width;
global m2_index_p;
global m2_id_width;

global G_MaxModels ;

for ap=1:nGroups    % ap means!

    % *******************
    %  validation.
    % *******************
    rowValid = 1;
    indexEnd   =  m3_header_width + m3_ids*m3_id_width;

    indexStart =  m3_header_width + m3_index_strategy2_y_estimate;
    YIdx       =  indexStart : m3_id_width : indexEnd; 
    YEstVector   =  m3(ap, YIdx);

    indexStart =  m3_header_width + m3_index_strategy2_y_estimate2;
    YIdx       =  indexStart : m3_id_width : indexEnd; 
    YModelVector  =  m3(ap, YIdx);

    if (sum(YEstVector>0) < m3_ids*0.7)
        rowValid = 0;
    end
    if (sum(YModelVector>0) < m3_ids*0.7)
        rowValid = 0;
    end

    % *******************
    %  temp.
    % *******************
    if (ap > 3)
        rowValid = 1;
    end

    %m3(ap, m3_index_global_valid) = rowValid;
    if (rowValid < 1)
        continue;
    end


    for iModel=1:G_MaxModels;
        indexStart =  m3_header_width + m3_index_strategy2_b0 + iModel ;  % y-est modelx (not sort)
        indexL = indexStart:m3_id_width:indexEnd ;
        rankList   =  m3(ap, indexL);

        % **************************************************
        %  special for model2
        %    
        % **************************************************
        if (iModel == 112)
            onlyPart1 = 0;
            onlyPart2 = 0;
            if (onlyPart1 == 1)
                rankList(m3_ids/2:m3_ids) = 0;
            end
            if (onlyPart2 == 1)
                rankList(1:m3_ids/2) = 0;
            end
        end


        % **************
        %  ....
        % **************
        bestNum = sum(rankList==1)


        % **************************************************
        %  special for model1
        %    random
        % **************************************************
        if (iModel == 1)
            rankList = rand(size(rankList));
        end
 

        n = length(rankList);
        [a I] = sort(rankList);
        tmp(I(:)) = 1:n;
        iModel
        sortRankList = tmp



        % **************************************************
        %  special for model3 
        %    set same zeros for sortRankList[]
        % **************************************************
        if (iModel == 3)
            %rankList   =  rand(1,m3_ids);
            % ******************************************
            %  set to 0 for the non-predictable ones.
            % ******************************************
            clearOnes = 0;
            errList = 0;
            for sid=1:m3_ids


break;

                %index_yReal    =  m3_header_width + (sid-1)*m3_id_width + m3_index_strategy2_y;
                index_yReal    =  m3_header_width + (sid-1)*m3_id_width + m3_index_strategy2_b0_pvalue - 2;
                index_yModel3  =  m3_header_width + (sid-1)*m3_id_width + m3_index_strategy2_b0_pvalue - 0;
                apStart = ap - 5;
                apEnd   = ap - 1;
                if (apStart < 1)
                    break;
                end
                data_yReal   = m3(apStart:apEnd, index_yReal);
                data_yModel3 = m3(apStart:apEnd, index_yModel3);
    
                % ******************************************
                %  set to 0 for the non-predictable ones.
                % ******************************************
                %data_yReal = data_yReal/m3_ids;
                data_yReal = data_yReal;
                errPredict = mean( abs(data_yReal - data_yModel3))
                errList(sid) = errPredict;
                if (0)
                %if (errPredict > 0.4)
                    data_yReal
                    data_yModel3
                    clearOnes++;
                    sortRankList(sid) = 0;
                end
            end
            ap
            clearOnes 
            errList
            errMean = mean(errList)
            sortRankList(errList > errMean) = 0;
            %sortRankList(errList > 0.20 ) = 0;
        end



        % ******************
        %   stream to M3
        % ******************
        indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue + iModel ;  % y-est modelx (sort)
        indexL = indexStart:m3_id_width:indexEnd ;
        m3(ap, indexL) = sortRankList ./ n ;



        % ******************
        %   boost 
        % ******************
        if (iModel==15)
            indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue + 4;
            indexL = indexStart:m3_id_width:indexEnd ;
            r4 = m3(ap, indexL);

            indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue + 5;
            indexL = indexStart:m3_id_width:indexEnd ;
            r5 = m3(ap,indexL);

            r5 = r5 + r4;
            %r5(r4<0.5) = 0;
            %r5 = r5 + r4;
            iiindex = r4>0.99;
            r5(iiindex) = r5(iiindex) + r4(iiindex);
            r5


            rankList = r5;
            tmp=0;
            [a I] = sort(rankList);
            tmp(I(:)) = 1:n;
            sortRankList = tmp;
            m3(ap,indexL) = sortRankList ./ n;
    
        end



        % ******************
        %   random 
        % ******************
        if (iModel==113) %specifyModel)
            indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue + 3;
            indexL = indexStart:m3_id_width:indexEnd ;
            r3 = m3(ap,indexL)

            indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue + 2;
            indexL = indexStart:m3_id_width:indexEnd ;
            rx = m3(ap,indexL);

            %r3 = rand(1,n);
            %r3 = r3+r2;
            %rx(r3>0.9) = 0;
            %rx(r3<0.4) = 0;
             removeIdx = rx<0.7;
             removeNums = sum(removeIdx)
             %rx(removeIdx) = 0;
            % r3(removeIdx) = 0;
             %r3 =rx + r3;
             %r3 = exp(rx.^2) + exp(r3.^2);

            if (1)
                r3(r3<0.99) = 0;
                rx(rx<0.99) = 0;
                rr = r3 + rx;
                %for ii=1:len
                %end
                r3 = rr
            end

            rankList = r3;
            tmp=0;
            [a I] = sort(rankList);
            tmp(I(:)) = 1:n;
            sortRankList = tmp;

            % write back
            indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue + 3;
            indexL = indexStart:m3_id_width:indexEnd ;
            m3(ap,indexL) = sortRankList ./ n;
        end
    end

end



% ********************************************
%  make the roi and real-y for each unit. 
% ********************************************
for ap=1:nGroups
    for sid=1:m3_ids
       %thisRoi  = M3_getRoiOfFuture(ap, sid);
       %thisRoi  = M3_getRoiOfFuture_dynSell(ap, sid);
       [thisRoi timeCost]  = M3_getRoiOfFuture_dynSell_avg(ap, sid);
       %[thisRoi timeCost]  = M3_getRoiOfFuture_dynSell_avg_afternoon(ap, sid);
       %[thisRoi timeCost]  = M3_getRoiOfFuture_avline(ap, sid);


       index_ = m3_header_width + m3_index_strategy1_roi        + (sid-1)*m3_id_width;
       m3(ap,index_)  = thisRoi;

       index_ = m3_header_width + m3_index_strategy2_indicator1 + (sid-1)*m3_id_width;
       m3(ap,index_)  = timeCost;
    end

    ap
    indexRoi_ap = m3_header_width + m3_index_strategy1_roi + ((1:m3_ids)-1)*m3_id_width;
    roiArray = m3(ap, indexRoi_ap)
    % roiArrayRegular = M3_regular(roiArray)

    rankList = roiArray;
    n = length(rankList);
    [a I] = sort(rankList);
    tmp(I(:)) = 1:n;
    roiArrayRegular = tmp/n
    
    indexYreal_ap = m3_header_width + m3_index_strategy1_y + ((1:m3_ids)-1)*m3_id_width;
    m3(ap, indexYreal_ap) = roiArrayRegular;
end


endfunction





######################################################## FILE.
function  [] = M3_writeBack()
######################################################## FILE.
global m3 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dlmwrite("data_m3.edat", m3, "\t", 'precision', '%2.9f');

endfunction



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ytest rosv] = svr_trainer2(xdata,ydata, xtest, C, epsilon, qpIter, kernel, varargin)
% SVR  Utilises Support Vector Regression to approximate 
%           the functional relationship from which the
%           the training data was generated.
%  Function call:
%
%    svrobj = svr_trainer(x_train,y_train,c,epsilon,kernel,varargin);
%    The training data, x_train and y_train must be column vectors.
%
%  Example usage:
%
%    svrobj = svr_trainer(x_train,y_train,400,0.000000025,'gaussian',0.5);
%    y = svrobj.predict(x_test);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(kernel,'gaussian')
        lambda = varargin{1};
        kernel_function = @(x,y) exp(-lambda*norm(x-y,2)^2);
    elseif strcmp(kernel,'spline')
        %kernel_function = @(a,b) prod(arrayfun(@(x,y) 1 + x*y+x*y*min(x,y)-(x+y)/2*min(x,y)^2+1/3*min(x,y)^3,a.feature,b.feature));
        kernel_function = @(a,b) prod(arrayfun(@(x,y) 1 + x*y+x*y*min(x,y)-(x+y)/2*min(x,y)^2+1/3*min(x,y)^3,a,b));
    elseif strcmp(kernel,'periodic')
        l = varargin{1};
        p = varargin{2};
        kernel_function = @(x,y) exp(-2*sin(pi*norm(x.feature-y.feature,2)/p)^2/l^2); 
    elseif strcmp(kernel,'tangent')
        a = varargin{1};
        c = varargin{2};
        kernel_function = @(x,y) prod(tanh(a*x.feature'*y.feature+c)); 
    elseif strcmp(kernel,'corr')
        kernel_function = @kernel_func_corr;
    elseif strcmp(kernel,'poly')
        kernel_function = @(x,y) (dot(x,y) +1)^2;
    end
    
    
    Options.MaxIter = qpIter
    ntrain = size(xdata,1);
    alpha0 = zeros(ntrain,1);
    
    % *********************************
    %  meta parameters.
    % 
    % *********************************
    %lambda = 0.5;
    %C = 100000;
    %epsilon = 0.005;


    % *********************************
    % Set up the Gram matrix for the 
    % training data.
    % *********************************
    for i=1:ntrain
    	for j=1:i
            xi= xdata(i,:);
            xj= xdata(j,:);
            K(i,j) = kernel_function(xi, xj);
            K(j,i) = K(i,j);
        end
    end
    K
    detK = det(K)
    if (detK ==0)
        svmEcept = -1
        ytest = 0;
        rosv = 1;
        return;        % Note!!!!!!!
    end



    
    % *********************************
    % Set up qp. 
    % *********************************
    H   = [K -K; -K K];
    f   = [ epsilon*ones(ntrain,1) - ydata; epsilon*ones(ntrain,1) + ydata];
    lb  = [ zeros(2*ntrain,1);];
    ub  = [C*ones(2*ntrain,1);];
    Aeq = [ones(ntrain,1); -1*ones(ntrain,1);]'; % row vector.
    beq = 0;

    
    % *********************************
    % Train the SVR by optimising the  
    % dual function ie. find a_i's 
    % *********************************
    
    %options = optimoptions('quadprog','Algorithm','interior-point-convex');
    %figure; imagesc(M); title('Inner product between training data (ie. K(x_i,x_j)'); xlabel('Training point #'); ylabel('Training point #');
    xx0 = zeros(2*ntrain,1);
    [z obj info] = qp(xx0, H,f,Aeq,beq,lb,ub, Options);
    iter=info.solveiter
    if (info.info == 0)
    else
        maxIterErr=info.info
        [ntest _nn] = size(xtest);
        %ytest =  zeros(1, ntest);
        ytest =  -7.777*ones(1, ntest);
        rosv = 1;
        return;        % Note!!!!!!!
    end
    
    alpha1 = z(1:ntrain);
    alpha2 = z(ntrain+1:2*ntrain);
    alpha  = [alpha1 alpha2 ]
    %figure; stem(alpha); title('Visualization of the trained SVR'); xlabel('Training point #'); ylabel('Weight (ie. alpha_i - alpha_i^*)');


    % *********************************
    % Calculate b  and rate of SV
    % *********************************
    nSV = 0;
    for i=1:ntrain
        b(i) = ydata(i) + epsilon;
        %xi = xdata(i, :);
        for j = 1:ntrain
            %xj = xdata(j, :);
            %b(i) = b(i) - (alpha1(j) - alpha2(j))*kernel_function(xi, xj);
            b(i) = b(i) - (alpha1(j) - alpha2(j))*K(i,j);
        end


        zeroNum = e^-8;
        if (and(alpha1(i)<zeroNum, alpha2(i)<zeroNum))
            %nInner++;
        else
            nSV++;
        end
    end
    b
    b = mean(b);
    rosv = nSV/ntrain
    
    
    % *********************************
    %  predict.
    % *********************************
    [ntest _nn] = size(xtest);
    for i=1:ntest
        ytest(i) = b;
        xx = xtest(i, :);
        for j=1:ntrain
           factor = alpha1(j) - alpha2(j);
           if (abs(factor) < e-10)
               continue;
           end

           xj = xdata(j,:);
           ytest(i) = ytest(i) + factor*kernel_function(xx, xj);
        end
    end

    return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




######################################################## FILE.
function  [] =  M3_postModel3()
######################################################## FILE.
global m3 ;
global m3_ids ;
global nGroups;
global m2_groupSize;

global m3_header_width;
global m3_index_strategy2_y_estimate;
global m3_index_strategy2_y_estimate2;
global m3_index_strategy1_roi;
global m3_id_width;
global m3_index_date1     ;    % header area
global m3_index_date2     ;    % header area
global m3_index_global_valid;  % header area
global m3_index_envStress ;%= 3;    % header area
global m3_index_r1        ;%= 5;   % header area
global m3_index_global_valid;
global m3_index_market    ;   % header area
global m3_index_market_1  ;   % market roi index
global m3_index_market_2  ;   % market roi index
global m3_index_market_3  ;   % header area

global m3_index_strategy2_b1;
global m3_index_strategy2_b0;
global m3_index_strategy2_b0_pvalue;
global m3_index_strategy2_y ;
global m3_index_strategy1_y ;



global m2;
global m2_header_width;
global m2_index_p;
global m2_index_v;
global m2_index_jbm_tv;
global m2_id_width;
            
global m2_header_p_index;   % header area

global stockIdStr;
global G_pg;
global G_ModelSelect_toFile;

global m3_index_strategy2_b7_pvalue  ;
global m3_index_strategy2_b8_pvalue  ;

global stockBriefIDStr;
global stockBriefIPO;
global stockBriefGBLT;
global stockBriefNetProfit;

for beta=1:1:1

    indexEnd   =  m3_header_width + m3_ids*m3_id_width;
    randList =  rand(1, m3_ids)
    if (1)
    for ap=1:nGroups


            %%%%%%%%%%%%%%%%%%%%
            %% m2 timeseries.
            %%%%%%%%%%%%%%%%%%%%
            thisAPidx  = (m2(:,2) == ap);
            timeSeries = m2(thisAPidx, 1); % set timestamp of each ap.
            tt = timeSeries(1);
            st = localtime(tt);
            ymd1 = (st.year+1900)*10000 + (st.mon+1)*100 + st.mday;
            m3(ap, m3_index_date1) = ymd1;
            tt = timeSeries(length(timeSeries));
            st = localtime(tt);
            ymd2 = (st.year+1900)*10000 + (st.mon+1)*100 + st.mday;
            %m3(ap, m3_index_date2) = ymd2 + ceil(st.wday)/10;
            m3(ap, m3_index_date2) = tt; %


            %%%%%%%%%%%%%%%%%%%%
            %% market.
            %%%%%%%%%%%%%%%%%%%%
            startidx  =  m2_header_width + m2_index_p - 2;
            endidx    =  m2_header_width + m3_ids*m2_id_width;
            idxPrawBlock   =  m2(thisAPidx, startidx:m2_id_width:endidx);
            idxPraw  = idxPrawBlock(size(idxPrawBlock,1), :);
            
            market = sum(idxPraw) / 10000;
            m3(ap, m3_index_market ) = market;

            lastAp = ap-1;
            if (lastAp>0)
                lastMarket = m3(lastAp, m3_index_market ) ;
                roiRatio = (market - lastMarket )/lastMarket;
                lastRoiRatio = m3(lastAp, m3_index_market_1 );
                m3(ap, m3_index_market_1 ) = roiRatio;
                m3(ap, m3_index_market_2 ) = mean([roiRatio lastRoiRatio]);
            end





        % **********************************
        %  generate feature of model3
        % **********************************
        thisGrIndex = (m2(:,2) == ap);
        PIndex = m2(thisGrIndex, m2_header_p_index)
        

        % **********************************
        %  get the style of market.
        % **********************************
        countPositive = 0;
        countValid    = 0;
        marketRoi    = 0;
        for sid=1:m3_ids
                m2IndexBase = m2_header_width + (sid-1)*m2_id_width;
                PSelfRaw  = m2(thisGrIndex, m2IndexBase + m2_index_p - 2);   % Note !!!!!!!!!!!!!
                n = length(PSelfRaw);

                if (PSelfRaw(n) > e^-9)
                    countValid++;
                else
                    continue;
                end

                %vsValue = PSelfRaw(n-24);
                vsValue = PSelfRaw(n-51);
                if (PSelfRaw(n) > vsValue )
                    countPositive ++;
                end

                if (vsValue > 0)
                    marketRoi  += (PSelfRaw(n) - vsValue)/vsValue;
                end
        end
        ratioPositive = countPositive/countValid;
        marketRoi = marketRoi/countValid;
        m3(ap, (m3_index_market_3)) = ratioPositive;
        printf ("ratioPositive %f, ap=%d validCount=%d marketRoi=%f\n", ratioPositive, ap, countValid, marketRoi);


        % **********************************
        %  init feature list.
        % **********************************
        featureList0_qd(ap, :) =  zeros(1, m3_ids);
        featureList1_qd(ap, :) =  zeros(1, m3_ids);
        featureList2_qd(ap, :) =  zeros(1, m3_ids);
        featureList3(ap, :) =  zeros(1, m3_ids);
        featureList4(ap, :) =  zeros(1, m3_ids);
        featureList5_hcl(ap, :) =  zeros(1, m3_ids);
        featureList6_x(ap, :) =  zeros(1, m3_ids);
        featureList7_cje(ap, :) =  zeros(1, m3_ids);
        featureList8_mcao_up(ap, :) =  zeros(1, m3_ids);
        featureList9_mcao_do(ap, :) =  zeros(1, m3_ids);
        featureList10_activity(ap, :) =  zeros(1, m3_ids);
        featureList11_cjeRal(ap, :)   =  zeros(1, m3_ids);
        featureList12_mrv(ap, :)   =  zeros(1, m3_ids);


        % **********************************
        %  compute Q indicator for each sid.
        % **********************************
        for sid=1:m3_ids
            m2IndexBase = m2_header_width + (sid-1)*m2_id_width;
            P = m2(thisGrIndex, m2IndexBase + m2_index_p - 1);   % Note !!!!!!!!!!!!!
            P = m2(thisGrIndex, m2IndexBase + m2_index_p );      % Note !!!!!!!!!!!!!
            V = m2(thisGrIndex, m2IndexBase + m2_index_v - 4);   % Note: magic number.
            n = length(P);


            % ***********************
            %   set zero(clean).
            % ***********************
            m3IndexBase = m3_header_width + (sid-1)*m3_id_width;
            m3(ap, m3IndexBase + m3_index_strategy2_b7_pvalue)  = 0;
            m3(ap, m3IndexBase + m3_index_strategy2_b8_pvalue)  = 0;
            m3(ap, m3IndexBase + m3_index_strategy2_b8_pvalue+1)= 0;
            m3(ap, m3IndexBase + m3_index_strategy2_b8_pvalue+2)= 0;


            % ************************
            %   space test policy.
            % ************************
            if (sid > m3_ids/2)
                %onlyPart1 = 1;
                %continue;
            end


            if (countValid > 1000)
            %if (randList(sid) > 1/18)
            if (randList(sid) > 0.5)
                %continue;
            end
            end


            % **************
            %   valid check.
            % **************
            if (sum(V==0) > n/10)
                %continue;
            end


            % ************************************
            %   compute average volume per 5 min.
            % ************************************
            VInt= m2(thisGrIndex, m2IndexBase + m2_index_v - 3);  % Note: magic number.
            %toV = VInt(n) - VInt(n-20);
            %avV = toV/20;
            toV = VInt(n);
            avV = toV/n;
            if (avV==0)
                continue;
            end



            % **************
            %   Q factor.
            % **************
            qIndicator1 = 0;
            qIndicator2 = 0; 
            qIndicator3 = 0; 
            qIndicator_A1 = 0;
            qIndicator_P2 = 0; 
            qIndicator_A2 = 0; 
            mcaoNum    = 0;
            if (1)
                PSelfRaw  = m2(thisGrIndex, m2IndexBase + m2_index_p - 2); % Note !!!!!!!!!!!!!
                %PSelf  = P;                                                % Note !!!!!!!!!!!!!
                PSelf        = m2(thisGrIndex, m2IndexBase + m2_index_p - 1); % Note !!!!!!!!!!!!!
                %PSelf        = PSelfRaw
                PSelfInvert  = -PSelf;

                PSelf_jbm_tv = m2(thisGrIndex, m2IndexBase + m2_index_jbm_tv);

                % *******************************
                %   jbm
                % *******************************
                jbm_tv = PSelf_jbm_tv(1);
                jbm_gb = PSelf_jbm_tv(1)/PSelfRaw(1);

                % *******************************
                %  same special case: price stay.
                % *******************************
                %if (sum(diff(PSelfRaw)==0) > n/2)
                %if (sum(diff(PSelfRaw(95:100))==0) > 4)
                %if (sum(diff(PSelfRaw(95:100))==0) > 3)
                %if (PSelfRaw(99) == PSelfRaw(100))
                %    tooMuchSameP = 1
                %    continue;
                %end


                % *******************************
                %  filter ...
                % *******************************
                idsStr = stockIdStr(sid);
                if (idsStr>300000 && idsStr<400000)
                    continue;
                end
                printf("sid:%d\t jbm_tv: %f, jbm_gb:%f\n", idsStr, jbm_tv, jbm_gb);


                % *******************************
                %  filter ...
                %  only consider cixin gu.
                % *******************************
                thisTs = m2(thisGrIndex, 1);
                %thisTs = thisTs(1);
                hsl = 0.;
                %if (thisTs - stockBriefIPO(sid) > 3600*24*365*3)
                    %continue;
                %end
                %if (stockBriefGBLT(sid) > 0)
                if (1)
                    thisV = VInt(n) - VInt(n-25);
                    hsl = (thisV*100)/stockBriefGBLT(sid);
                    %hsl = (toV*100)/stockBriefGBLT(sid);
                    %if (hsl < 0.01)
                    %    %continue;
                    %end
                end
                printf("sid:%d,\t thisTs:%d,\t IPO:%d,\t hsl:%f\n", idsStr, thisTs, stockBriefIPO(sid), hsl);



                % *******************************
                %   consider jbm.
                %    discard it in case of negative profit 
                % *******************************
                if (stockBriefNetProfit(sid) < 10)
                    continue;
                end


                % *******************************
                %   market value management.
                % *******************************
                cje = toV * 10^6 * PSelfRaw(1)
                %if (cje > 10^8*20)
                if (cje > 10^8*10)
                    cjeTooBig = 1
                    %continue;
                end
                if (cje < 10^8*0.5)
                    cjeTooSmall = 1
                    %continue;
                end


                gagerCount_neg = 0;
                positiveVol = 0;
                qValueList = zeros(1,n);
                qValueList2= zeros(1,n);
                qValueList3= zeros(1,n); % negative list1
                qValueList4= zeros(1,n);
                qValueList5_upv= zeros(1,n);
                qValueList6_dov= zeros(1,n);
                countUPV=0;
                countDOV=0;
    
                %for i=80:n   % Note: the range is sensetive for result. 2:n
                for i=4:n   % Note: the range is sensetive for result. 2:n
                    gapIndex      = PIndex(i) - PIndex(i-1) ;
                    gapIndex_pre  = PIndex(i-1) - PIndex(i-2) ;
                    gapIndex_pre2 = PIndex(i-2) - PIndex(i-3) ;
                    gapIndividual = PSelf(i)  - PSelf(i-1) ;
                    gapIndividual_pre = PSelf(i-1)  - PSelf(i-2) ;
                    %gapIndividual_pre2= PSelf(i-2)  - PSelf(i-3) ;
                    %gapIndividual_enhance     = PSelf(i)  - PSelf(i-1) ;
                    %gapIndividual_enhance_pre = PSelf(i-1)- PSelf(i-2) ;
                    gapIndividual = gapIndividual * PSelfRaw(1) / PSelfRaw(i-1);


                    %pCoef = 50*1;
                    %pCoef = 50*2.0*(2-ratioPositive);
                    %pCoef = 50*3.0;
                    pCoef = 50*2.0;
                    vCoef = -0.05;
                    vCoef = 0.0;
                    vCoef = -0.01;
                    overCoef = 2;
                    overCoef = 3.14159;
                    %eCoef = 0.0005
                    eCoef = 0.001
                    fangzhi = 0.0005;
                    fangzhi = 0.0000;


                    if (1) 
                        vol = VInt(i) - VInt(i-1);
                        if (vol < 0)
                            findVolException = 1
                            vol
                            vol = 0;
                        end
                        %if (vol/avV > 3)
                        %    mcaoNum ++;
                        %end
                        if ( 1*(gapIndividual) > 0.007)
                            %mcaoNum ++;
                            qValueList2(i) = abs(gapIndividual);
                            %qValueList2(i) = 1;
                        end
                        if (-1*(gapIndividual) > 0.007)
                            %mcaoNum ++;
                            qValueList3(i) = abs(gapIndividual);
                            %qValueList3(i) = 1;
                        end


                        if ((i>75 || (i>25&&i<75)) && gapIndividual > 0.003)
                            countUPV ++; 
                            qValueList5_upv(i)   = vol / gapIndividual ;
                        end
                        if ((i>75 || (i>25&&i<75)) && gapIndividual < -0.003)
                            countDOV ++;
                            qValueList6_dov(i)   = -vol / gapIndividual;
                        end
                    end

                    
                    %if (i>75) 
                    %    qValueList(i) = PSelf(i) - PSelf(74);
                    %end
                    %baseValue = PSelf(74) - 0.2;
                    baseValue =  0 - 0.2;

                    if (1)
                    if (gapIndividual>fangzhi)
                        pFactor = (gapIndividual)*pCoef
                        vFactor = vol/(avV);
                        vFactor = min(7,vFactor);
                        vFactor2=vFactor;
                        vFactor=eCoef*exp(vCoef*vFactor);

                        % same principle with avV, but avV will be filter by a max-value.
                        positiveVol += vol;

                        qValue = (1+pFactor)*vFactor;
                        %qValue = PSelf(i) - baseValue;
                        qValueList(i)   += qValue;
                        qValueList4(i)  += pFactor*vFactor;
                    end
                    if (gapIndividual<-fangzhi)
                        pFactor = -1 * (gapIndividual)*pCoef;
                        vFactor = vol/(avV);
                        vFactor = min(7,vFactor);
                        vFactor=eCoef*exp(vCoef*vFactor);

                        qValue = (-1)*(1+pFactor)*vFactor;
                        %qValue = -1*(PSelf(i) - baseValue);
                        qValueList(i)   += qValue;
                        %qValueList4(i)  += -1;
                    end


                    % ********************
                    %  overlap 
                    % ********************
                    if (gapIndividual>fangzhi && gapIndividual_pre>fangzhi)
                        %pFactor = (gapIndividual+gapIndividual_pre)*pCoef
                        pFactor = (gapIndividual)*pCoef
                        vol2 = VInt(i) - VInt(i-2);
                        vFactor = vol2/(2*avV);
                        vFactor = min(7,vFactor);
                        vFactor2=vFactor;
                        vFactor=eCoef*exp(vCoef*vFactor);

                        qValueList(i)  += (1+pFactor)*vFactor*overCoef;
                        %qValueList(i)  += (PSelf(i) - baseValue);
                        qValueList4(i) += 1;
                        %qValueList4(i)  += 2*pFactor*vFactor;
                    end

                    if (gapIndividual<-fangzhi && gapIndividual_pre<-fangzhi)
                        pFactor = -1*(gapIndividual+gapIndividual_pre)*pCoef
                        %pFactor = -1*(gapIndividual)*pCoef
                        vol2 = VInt(i) - VInt(i-2);
                        vFactor = vol2/(2*avV);
                        vFactor = min(7,vFactor);
                        vFactor=eCoef*exp(vCoef*vFactor);

                        %qValueList(i)  += -2*(1+pFactor)*vFactor*overCoef;
                        qValueList(i)  += -1*(1+pFactor)*vFactor*overCoef;
                        %qValueList(i)  += -1*(PSelf(i) - baseValue);
                        %signal2Count = signal2Count + 2.3;
                    end
                    end


                    if (gapIndividual >0)
                        if (i>76 &&  gapIndividual>0.04)
                            %thisTs = m2(thisGrIndex, 1);
                            printf("find mengla:(ts=%d,\t i=%d,\t sstr=%d) ap=%d, sid=%d\n", thisTs(i), i, idsStr, ap, sid);
                        end
                    end
                end
    
                
                printf("ap=%d, sid=%d, s1,s2,s3 -> (%d, %d, %d)\n", ap, sid, 0,0,0);
                %fftP = abs(fft(PSelf));
                %[a_ b_] = max(fftP(2:n/2));
                %fft_qInd1 = a_ 
                %fft_qInd2 = b_;
                %fft_qInd3 =  std(PSelf)
            end


            % ----------------------------------------------
            cjeAM = VInt(n) - VInt(n-20);
            cjeBefore = VInt(n-20) - VInt(1);
            cjeRal = 0;
            if (ap>1 && cje>10^8*1)
                %if (tmp>0)
                if (cjeBefore>0)
                    cjeRal = cjeAM/cjeBefore;
                end
            end


            % ----------------------------------------------
            avl = M3_getAver(PSelf, 5);
            activity = M3_getActivity(PSelf, avl);


            % ----------------------------------------------
            bounce = 0
            stage2roi = 0;
            if (PSelfRaw(1)>0 && PSelfRaw(76)>0)
                stage1roi = (PSelfRaw(75)  - PSelfRaw(1)) /PSelfRaw(1)
                stage2roi = (PSelfRaw(100) - PSelfRaw(76))/PSelfRaw(76)
                if (stage1roi>0.1 && stage2roi<0)
                    bounce = exp(-stage2roi)
                end
            end
            stage2roi


            % ----------------------------------------------
            newqd = 0;
            qd3= 0;
            qd4= 0;
            newqd2= 0;
            stdvar = 0;
            nearstRoi = 0;
            throughRoi = 0;
            lastdayZhangfu = 0;
            jumpRoi = 0;
            %startPos = 86
            %startPos = 82
            %lastdayPos =49;
            %lastdayPos =89;
            %startPos = 50;
            %startPos = 88;
            startPos = 90;
            endPos = startPos+10;
            lastdayPos =74;
            lastlastdayPos = lastdayPos - 50;
            if (PSelfRaw(startPos)>0)
                nearstRoi = (PSelfRaw(endPos) - PSelfRaw(startPos))/PSelfRaw(startPos)
                stdvar = std(PSelfRaw(startPos:n))/PSelfRaw(n)
            end
            if (PSelfRaw(lastdayPos)>0)
                throughRoi = (PSelfRaw(endPos) - PSelfRaw(lastdayPos))/PSelfRaw(lastdayPos)
                jumpRoi = (PSelfRaw(lastdayPos+2) - PSelfRaw(lastdayPos))/PSelfRaw(lastdayPos)
                if (PSelfRaw(lastlastdayPos) > 0)
                    lastdayZhangfu = (PSelfRaw(lastdayPos) - PSelfRaw(lastlastdayPos)) / PSelfRaw(lastlastdayPos);
                end
            end

            %MV = VInt(startPos) - VInt(lastdayPos);
            MV = VInt(endPos) - VInt(startPos);
            %MV_before = VInt(100-50) - VInt(startPos-50);
            %MV_before = VInt(100) - 0;
            %MV_before = VInt(lastdayPos+2) - VInt(lastdayPos);
            MV_before = VInt(endPos) - VInt(lastdayPos);
            MRV = 0;
            if (MV>0 && MV_before>0)
                MRV = MV/MV_before
            end


            % uplimit == 0.02 will be worse.
            %if (nearstRoi<0.01 && nearstRoi>=0.00)
            %if (nearstRoi<0.01 && nearstRoi>=-0.01)
            %if (nearstRoi<0.01 && nearstRoi>= -0.01)
            %if (nearstRoi<0.01)
            %if (throughRoi<0.01)
            %if (throughRoi<0.05 && nearstRoi>-0.01 && nearstRoi<0.01)
            %if (throughRoi<0.05 && nearstRoi>0.0)
            %if (throughRoi<0.05)
            %if (throughRoi<0.03 && nearstRoi>0.0)
            %if (throughRoi<0.07)
            if (jumpRoi>-0.02 && throughRoi<0.05 && countDOV>1 && countUPV>1)
            %if (throughRoi<0.05 && nearstRoi>0.0)
            %if (throughRoi<0.07 && nearstRoi>0.0)
            %if (throughRoi<0.08 && nearstRoi>0.0 && nearstRoi<0.02)
                %stepRoir= (PSelfRaw(startPos) - PSelfRaw(60))/PSelfRaw(60)
                %newqd = sum(qValueList(startPos:n ))
                %newqd = max(0, sum(qValueList(startPos:n )))
                newqd =  0 + sum(qValueList(startPos:endPos))
                %newqd2=      sum(qValueList4(startPos:endPos))
                % stdvar = std(PSelf(startPos:n))
                % stdvar = PSelf(n) - min(PSelf(startPos:n))
                %stdvar = M3_getHcl(PSelfRaw(startPos:n))
                %stdvar = std(PSelfRaw(startPos:n))/PSelfRaw(n)

                if (newqd<0)
                    newqd=0;
                end

                meanv = mean(PSelfRaw(startPos:endPos));
                if (meanv>0)
                    vv = (PSelfRaw(endPos) - meanv)/meanv;
                    if (vv<0)
                        newqd = 0;
                        newqd2= 0;
                        %newqd = -vv;
                    else
                        newqd = newqd;
                        newqd2= vv;
                        %newqd = 0;
                    end
                end


                %if (countDOV>0 && countUPV>0)
                %if (1)
                if (newqd>0)
                    unitDOV = sum(qValueList6_dov) / countDOV;
                    unitUPV = sum(qValueList5_upv) / countUPV;
                    %qd3 = unitDOV / unitUPV;
                    %qd4 = 1/ qd3;
                    %printf ("[%2d,%6d] qd3=%f \t (%d,\t%d)unitDOV,unitUPV=%f,\t%f\n",ap, idsStr,  qd3, countDOV, countUPV, unitDOV, unitUPV);

                    %qthis = qValueList6_dov(75:n);
                    %qbefo = qValueList6_dov(25:75);
                    qthis = qValueList5_upv(75:n);
                    qbefo = qValueList5_upv(25:75);
                    qthisn = sum(qthis>0);
                    qbefon = sum(qbefo>0);
                    idxsel = qthis > 0;
                    qthisdov = qthis(idxsel)
                    if (qthisn>0 && qbefon>0)
                        %thisdayUnitDOV = sum(qthis) / qthisn;
                        %befodayUnitDOV = sum(qbefo) / qbefon;
                        thisdayUnitUPV = sum(qthis) / qthisn;
                        befodayUnitUPV = sum(qbefo) / qbefon;
                       
                        %qd3 = thisdayUnitDOV/befodayUnitDOV;
                        qd3 = thisdayUnitUPV/befodayUnitUPV;
                        qd4 = 1/ qd3;
                        printf ("[%2d,%6d] qd3=%f \t (%d,\t%d)unitDOV,unitUPV=%f,\t%f\n",ap, idsStr,  qd3, countDOV, countUPV, unitDOV, unitUPV);
                    else
                        newqd = 0;
                    end
                end
            end



            % ----------------------------------------------
            hsl2 = 0;
            %meanv = mean(PSelfRaw(startPos:endPos));
            %if (throughRoi<0.08 && meanv>0 && lastdayZhangfu>0.098)
            %if (throughRoi<0.08 && meanv>0 && lastdayZhangfu<-0.05)
            %if (throughRoi<0.08 && meanv>0 && lastdayZhangfu>0.04)
            %if (throughRoi<0.08 && meanv>0 && lastdayZhangfu<-0.03)
            %if (throughRoi<0.08 && meanv>0 && lastdayZhangfu>0.03)
            %if (throughRoi<0.08 && PSelfRaw(n)>meanv)
            if (throughRoi<0.08)
                hsl2 = hsl;
            end

  
            % ----------------------------------------------
            nearstRoi2 = 0;
            stdvar2 = 0;
            if (PSelfRaw(74)>0)
                nearstRoi2 = (PSelfRaw(100) - PSelfRaw(74))/PSelfRaw(74)
            end
            if (nearstRoi2<0.01)
                stdvar2 =  n + sum(qValueList(74:n ))
            end



            % ----------------------------------------------
            %hcl = -1*M3_getHcl(PSelfRaw(51:100));
            hcl = -1*M3_getHcl(PSelfRaw(startPos:endPos));
            hcl2= -1*M3_getHcl(PSelfRaw(26:startPos));



            % ----------------------------------------------
            featureList0_qd(ap, sid)  =  hcl2;
            %featureList0_qd(ap, sid)  =  sum(qValueList(80:n ))/n;
            %featureList1_qd(ap, sid)  =  sum(qValueList)/n;
            %featureList1_qd(ap, sid)  =  sum(qValueList4(80:n))/n;
            featureList1_qd(ap, sid)  =  qd4;
            %featureList2_qd(ap, sid)  =  signal1Count / signal2Count;
            %featureList2_qd(ap, sid)  =  signal1Count - signal2Count;
            %featureList2_qd(ap, sid)  =  sum(qValueList4(80:n ))/n;
            %featureList2_qd(ap, sid)  =  sum(qValueList(80:n )) * exp(-1000*(stage2roi-0.02)^2) / n;
            featureList2_qd(ap, sid)  =  newqd/1;
            %featureList3(ap, sid)     =  std(PSelf);
            %featureList3(ap, sid)     =  stdvar;
            featureList3(ap, sid)     =  qd3;
            featureList4(ap, sid)     =  PSelf(n) - PSelf(n/2);
            featureList5_hcl(ap, sid) =  hcl;
            %featureList5_hcl(ap, sid) =  hsl;
            %featureList5_hcl(ap, sid) =  M3_getHcl(PSelfRaw(lastdayPos:endPos));
            %featureList5_hcl(ap, sid) =  M3_getHcl(PSelfRaw);
            %featureList5_hcl(ap, sid) =  bounce;
            %featureList5_hcl(ap, sid) =  stdvar2;
            %featureList6_x(ap, sid)  =  jumpRoi;
            %featureList6_x(ap, sid)  =  newqd2/n;
            featureList6_x(ap, sid)  =  hsl2;
            featureList7_cje(ap, sid)       = cje;
            %featureList7_cje(ap, sid)       = jbm_tv;
            %featureList7_cje(ap, sid)       = jbm_gb;
            featureList8_mcao_up(ap, sid)    = sum(qValueList2(50:82))/n;
            featureList9_mcao_do(ap, sid)    = sum(qValueList3(50:82))/n;
            featureList10_activity(ap, sid) =  activity;
            featureList11_cjeRal(ap, sid)   =  cjeRal;
            featureList12_mrv(ap, sid)   =  newqd2;
            sid
            activity
            mcaoNum
            % ----------------------------------------------

        end %endof sid



        % **********************************
        %   the next will use (ap-3)
        % **********************************
        if (ap<=4)
            continue;
        end


        % **********************************************
        %   gen qList 
        %  -------------
        % 
        % **********************************************
        %featureTmp1 = featureList4(ap-0, :)    + featureList4(ap-1, :)    + featureList4(ap-2, :)    + featureList4(ap-3, :) ;  
        %featureTmp1 = featureList4(ap-0, :)    + featureList4(ap-1, :)    + featureList4(ap-2, :)   ;  
        featureTmp1 = featureList4(ap-0, :)     + featureList4(ap-1, :)     + featureList4(ap-2, :)     + featureList4(ap-3, :)  + featureList4(ap-4, :) ;  
        featureTmp2 = featureList7_cje(ap-0, :) + featureList7_cje(ap-1, :) + featureList7_cje(ap-2, :) + featureList7_cje(ap-3, :);
        featureTmp3 = featureList8_mcao_up(ap-0, :) + featureList8_mcao_up(ap-1, :) + featureList8_mcao_up(ap-2, :);
        featureTmp4 = featureList9_mcao_do(ap-0, :) + featureList9_mcao_do(ap-1, :) + featureList9_mcao_do(ap-2, :);
        featureTmp5 = featureTmp3 + featureTmp4;
        %ft1 = M3_regular( featureList7_cje(ap,:) );
        %ft2 = M3_regular( featureTmp1 );
        %featureTmp4 = M3_regular( ft1 .* ft2 )
        %numHuges = max(10, floor(ratioPositive*80))
        %numHuges_invert = max(10, 80 - numHots);
        %numHuges = 25
        numHuges = 200
        %numHuges = 50
        numHots = max(5, floor(ratioPositive*30))
        numHots_invert = max(5, 30 - numHots);
        %---------------------
        %qList2 = featureList0_qd(ap, :);
        qList2 = featureList2_qd(ap, :);
        q2Num = sum(qList2 >0)
        if (q2Num < 3)
            randq2 = 1
            qList2 = rand(size(qList2));
        end


        
        %---------------------
        %q3_fac1 = M3_regular( M3_rectangle(featureList0_qd(ap,:), 1, 3, 1) );
        %q3_fac2 = M3_regular( featureList5_hcl(ap,:) );
        %q3_fac2 = 1 - M3_regular( featureList3(ap,:) );
        %q3_fac2 = 1 - M3_regular( featureTmp1 );
        %qList3 = M3_addZeroZ([q3_fac1; q3_fac2]);
        %q3Num = sum(qList3 >0)
        %---------------------
        %q3_fac1 = M3_regular( M3_rectangle(featureList0_qd(ap,:), 1, 20, 1) );
        %q3_fac1 = M3_regular( M3_rectangle(featureList6_qdInvert(ap,:), 1, 3, 1) );
        %q3_fac1 = M3_regular( M3_rectangle(featureList2_qd(ap,:), 1, 3, 1) );
        %q3_fac2 = M3_regular( featureList9_mcao(ap,:) );
        %qList3 = M3_addZeroZ([q3_fac1; q3_fac2]);
        %--------------------- Hot one
        %qList3_fac1 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, numHuges,  1) );
        %qList3_fac2 = M3_regular( featureTmp1 );
        %--------------------- Hot one
        %qList3_fac1 = M3_regular( M3_rectangle(featureTmp1, 1, numHuges,  1) );
        %qList3_fac2 = M3_regular( featureList7_cje(ap,:) );
        %qList3 = M3_addZero(qList3_fac1, qList3_fac2);
        %q3Num = sum(qList3 >0)
        %--------------------- Hot one
        %qList3_fac1 = M3_rectangle_mannul(featureList2_qd(ap,:));
        %qList3_fac2 = M3_regular( featureList7_cje(ap,:) );
        %qList3_1 = M3_addZero(qList3_fac1, qList3_fac2);
        %q3Num_ = sum(qList3_1 >0)

        %qList3_fac1 = M3_regular( M3_rectangle(qList3_1, 1, 100,  1) );
        %qList3_fac2 = M3_regular( featureTmp1 );
        %qList3 = M3_addZero(qList3_fac1, qList3_fac2);
        %q3Num = sum(qList3 >0)
        %--------------------- Hot one
        %qList3_fac0 = M3_rectangle_mannul(featureList2_qd(ap,:));
        %qList3_fac1 = M3_regular( M3_rectangle(featureTmp1, 1, 100,  1) );
        %qList3_fac2 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, 100,  1) );
        %qList3_fac3 = M3_regular( featureList7_cje(ap,:) );
        %qList3 = M3_addZeroZ([qList3_fac0; qList3_fac1; qList3_fac2; qList3_fac3]);
        %q3Num = sum(qList3 >0)
        %--------------------- Hot one
        %qList3_fac0 = M3_rectangle_mannul(featureList2_qd(ap,:));
        %qList3_fac0 = M3_rectangle_mannul(featureList0_qd(ap,:));
        %qList3_fac1 = M3_regular( M3_rectangle(featureTmp1, 1, 300,  1) );
        %qList3_fac2 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, 1500,  1) );
        %qList3_fac3 = M3_regular( featureList2_qd(ap,:) );
        %qList3_fac3 = M3_regular( featureTmp1 );
        %qList3_fac3 = M3_regular( featureList2_qd(ap,:) );
        %qList3_fac3 = M3_regular( featureList2_qd(ap,:) );
        %qList3 = M3_addZeroZ([  qList3_fac2; qList3_fac3]);
        %q3Num = sum(qList3 >0)

        % compute avg-cje-100
        %---------------------
        q_cur_cje   = featureList7_cje(ap,:);
        %q_tmp       = M3_rectangle(featureTmp1, 1, 50,  1);
        q_tmp       = M3_rectangle(featureList2_qd(ap,:), 1, 10,  1);
        avg_cje_100 = sum(q_cur_cje(q_tmp>0)) / 10.;
        printf ("ap=%d avg_shizhi_100=%f\n", ap, avg_cje_100);
        %q_cje_ok = q_cur_cje(q_cur_cje>avg_cje_100*0.7 & q_cur_cje<avg_cje_100*1.3);

        %---------------------
        %LIMIT = 0.99;
        %q3_fac1 = q_cur_cje>avg_cje_100*(1-LIMIT);
        %q3_fac2 = q_cur_cje<avg_cje_100*(1+LIMIT);
        %q3_fac3 = M3_regular( featureList2_qd(ap,:) );
        %qList3 = M3_addZeroZ([q3_fac1; q3_fac2; q3_fac3]);
        %q3Num = sum(qList3 >0)
        %---------------------
        q3Num_up = sum(featureTmp3>0)
        q3Num_do = sum(featureTmp4>0)
        %q3_fac1 = M3_rectangle(featureTmp1, 0, 300,  0) ;
        %q3_fac2 = M3_rectangle(featureTmp1, 1, 300,  0) ;
        %q3_fac4 = M3_rectangle(featureTmp3, 0, 300, 1);
        %q3_fac5 = M3_regular(featureList2_qd(ap,:));
        %qList3 = M3_addZeroZ([q3_fac1; q3_fac2;  q3_fac4;  q3_fac5]);
        %q3Num = sum(qList3 >0)
        q3_fac1 = M3_rectangle(qList2, 1, 3,  1);
        q3_fac2 = M3_regular( featureList3(ap,:) );
        qList3 = M3_addZero(q3_fac1, q3_fac2);
        q3Num = sum(qList3 >0)

        %qList4 = featureList4(ap,:);
        %q4Num = sum(qList4 >0)
        %---------------------
        %qList4_= abs(featureList4(ap-0, :)) + ... 
        %         abs(featureList4(ap-1, :)) + ...
        %         abs(featureList4(ap-2, :)) + ...
        %         abs(featureList4(ap-3, :)) ;
        %---------------------
        %q4_fac1 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, floor(countValid*0.02),  1) );
        %q4_fac2 = M3_regular( featureTmp1 );
        %q4_fac3 = M3_regular( featureTmp2 );
        %qList4 = M3_addZeroZ([q4_fac1; q4_fac2; (1-ratioPositive)*q4_fac3])
        %q4Num = sum(qList4 >0)
        %---------------------
        %qList4_fac2 = M3_regular( featureList0_qd(ap,:) );
        %qList4 = M3_addZeroZ([qList4_fac2;]);
        %idxx = featureList7_cje(ap,:) < 5*10^8;
        %qList4(idxx) = 0;
        %---------------------
        %q4_fac1 = M3_regular( M3_rectangle(featureList0_qd(ap,:), 1, 3, 1) );
        %q4_fac2 = M3_regular( featureList9_mcao(ap,:) );
        %---------------------
        %q4_fac1 = M3_regular( M3_rectangle(featureTmp1, 1, 30, 1) );
        %q4_fac2 = M3_regular( featureList2_qd(ap,:));
        %qList4 = M3_addZeroZ([q4_fac1; q4_fac2]);
        %q4Num = sum(qList4 >0)
        %---------------------
        %qList4_fac1 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, 20,  1) );
        %qList4_fac2 = M3_regular( featureTmp1 );
        %qList4 = M3_addZero(qList4_fac1, qList4_fac2)
        %q4Num = sum(qList4 >0)
        %---------------------
        %qList4 = featureList12_mrv(ap, :);
        %q4Num = sum(qList4 >0)
        %---------------------
        %q4_fac1 = M3_regular( M3_rectangle(qList3, 1, 15,  1) );
        %q4_fac2 = M3_regular( featureList12_mrv(ap,:) );
        %qList4 = M3_addZero(q4_fac1, q4_fac2)
        %q4Num = sum(qList4 >0)
        %qList4 = featureList2_qd(ap, :);
        %qList4 = featureList3(ap, :);
        %---------------------
        %qList4_fac1 = M3_regular( M3_rectangle(featureTmp1, 1, numHuges,  1) );
        %qList4_fac2 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, 30,  1) );
        %qList4_fac3 = M3_regular( featureList0_qd(ap,:) );
        %qList4 = M3_addZeroZ([ qList4_fac2; qList4_fac3]);
        %q4Num = sum(qList4 >0)
        %---------------------
        %qList4 = M3_regular( featureList11_cjeRal(ap,:) );
        %qList4_fac1 = M3_regular( M3_rectangle(featureTmp1, 1, 700,  1) );
        %qList4_fac3 = M3_regular( featureList2_qd(ap,:) );
        %qList4 = M3_addZeroZ([ qList4_fac1; qList4_fac3]);
        %q4Num = sum(qList4 >0)
        %---------------------
        %qList4 = featureList1_qd(ap, :);
        %q4Num = sum(qList4 >0)
        %---------------------
        %q4_fac1 = M3_rectangle(featureTmp1, 1, 200,  1) ;
        %q4_fac2 = M3_rectangle(featureTmp2, 1, 200,  1) ;
        %q4_fac3 = M3_regular(qList2);
        %qList4 = M3_addZeroZ([q4_fac1; q4_fac2; q4_fac3]);
        %qList4 = featureList3(ap,:);
        %---------------------
        q4_fac1 = M3_rectangle(featureTmp5, 0, 300, 0);
        q4_fac2 = M3_rectangle(featureTmp5, 1, 300, 0);
        q4_fac3 = M3_regular(featureList2_qd(ap,:));
        qList4 = M3_addZeroZ([q4_fac1; q4_fac2; q4_fac3]);
        q4Num = sum(qList4 >0)



        %---------------------
        %qList5 = featureList0_qd(ap,:);
        %q5Num = sum(qList5 >0)
        %---------------------
        %qList5_fac1 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, numHots,  1) );
        %qList5_fac2 = M3_regular( featureTmp1 );
        %qList5 = M3_addZero(qList5_fac1, qList5_fac2)
        %---------------------
        %qList5 = qList3
        %q5Num = sum(qList5 >0)
        %---------------------
        %q5_fac1 = M3_regular( M3_rectangle(featureList12_mrv(ap,:), 1, 50,  1) );
        %q5_fac2 = M3_regular( featureList2_qd(ap,:) );
        %qList5 = M3_addZero(q5_fac1, q5_fac2);
        %---------------------
        %qList5 = featureList5_hcl(ap, :);
        %---------------------
        %q5_fac1 = M3_regular( M3_rectangle(qList3, 1, 6,  1) );
        %q5_fac2 = M3_regular( featureList10_activity(ap,:) );
        %qList5 = M3_addZero(q5_fac1, q5_fac2);
        %q5Num = sum(qList5 >0)
        %---------------------
        %qList5_fac0 = M3_regular( M3_rectangle(featureTmp1, 1, 180,  1) );
        %qList5_fac1 = M3_regular( M3_rectangle(featureList2_qd(ap,:), 1, 1800,  1) );
        %qList5_fac2 = M3_regular( featureList11_cjeRal(ap,:) );
        %qList5 = M3_addZeroZ([ qList5_fac0; qList5_fac1; qList5_fac2]);
        %q5Num = sum(qList5 >0)
        %---------------------
        %qList5_fac0 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, 500,  1) );
        %qList5_fac1 = M3_regular( M3_rectangle(featureTmp1, 1, 100,  1) );
        %qList5_fac1 = M3_regular( M3_rectangle(featureList6_x(ap,:), 1, 100,  1) );
        %qList5_fac1 = M3_regular( featureList6_x(ap,:) );
        %qList5_fac3 = M3_regular( featureTmp1 );
        %qList5 = M3_addZeroZ([ qList5_fac0; qList5_fac1; qList5_fac3]);
        %q5Num = sum(qList5 >0)

        %---------------------
        q_cje_normal  = exp (-1 .* ((q_cur_cje - avg_cje_100) / avg_cje_100 ).^2)
        %qList5 = M3_regular( featureList2_qd(ap,:) ) .* q_cje_normal;
        %q5Num  = sum(qList5 >0)
        %---------------------
        %LIMIT = 0.80;
        % q5_fac1 = q_cur_cje>avg_cje_100*(1-LIMIT);
        % q5_fac2 = q_cur_cje<avg_cje_100*(1+LIMIT);
        % q5_fac3 = M3_regular( featureList2_qd(ap,:) );
        % qList5 = M3_addZeroZ([q5_fac1; q5_fac2; q5_fac3]);
        %q5Num = sum(qList5 >0)
        %---------------------
        %q5_fac1 = M3_rectangle(featureTmp1, 1, 100,  1) ;
        %q5_fac2 = M3_regular( featureList3(ap,:) );
        %qList5 = M3_addZeroZ([ q5_fac1; q5_fac2]);
        %q5Num = sum(qList5 >0)
        %---------------------
        %qList5 = M3_regular(featureList5_hcl(ap,:));
        %q5Num = sum(qList5 >0)
        %q5_fac1 = M3_rectangle(qList2, 1, 6,  1);
        %q5_fac2 = 1-M3_regular(featureList5_hcl(ap,:));
        %qList5 = M3_addZeroZ([ q5_fac1; q5_fac2]);
        %q5Num = sum(qList5 >0)
        %---------------------
        %q5_fac1 = M3_rectangle(featureTmp1, 0, 200,  0) ;
        %q5_fac2 = M3_rectangle(featureTmp1, 1, 200,  0) ;
        %q5_fac5 = M3_regular(featureTmp4);
        %qList5 = M3_addZeroZ([q5_fac1; q5_fac2; q5_fac5]);
        %q5Num = sum(qList5 >0)
        q5_fac1 = M3_rectangle(featureTmp5, 0, 500, 0);
        q5_fac2 = M3_rectangle(featureTmp5, 1, 500, 0);
        q5_fac3 = M3_regular(featureList2_qd(ap,:));
        qList5 = M3_addZeroZ([q5_fac1; q5_fac2; q5_fac3]);
        q5Num = sum(qList5 >0)







        %qList6_fac1 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, floor(countValid*0.02),  1) );
        %qList6_fac2 = M3_regular( featureTmp1 );
        %qList6 = M3_addZero(qList6_fac1, qList6_fac2)
        %----------------------------------
        %q6_fac1 = M3_regular( M3_rectangle(featureList0_qd(ap,:), 1, 5,  1) );
        %q6_fac2 = M3_regular( featureList7_cje(ap,:) )
        %qList6 = M3_addZero(q6_fac1, q6_fac2);
        %----------------------------------
        %qList6 = featureTmp4;
        %q6Num = sum(qList6 >0)
        %----------------------------------
        %q6_fac1 = M3_rectangle(featureTmp3, 1, 100, 1);
        %q6_fac2 = 1  - M3_regular( featureList10_activity(ap,:) );
        %qList6 = M3_addZeroZ([q6_fac1; q6_fac2])
        %----------------------------------
        %q6_fac1 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, 1500,  1) );
        %q6_fac2 = M3_regular(featureList0_qd(ap,:));
        %qList6 = M3_addZeroZ([q6_fac1; q6_fac2])
        %q6Num = sum(qList6 >0)
        %----------------------------------
        %q6_fac1 = M3_regular( M3_rectangle(featureList0_qd(ap,:), 1, 3, 1) );
        %q6_fac2 = 1 - M3_regular( featureTmp1 );
        %qList6 = M3_addZeroZ([q6_fac1; q6_fac2]);
        %q6Num = sum(qList6 >0)
        %----------------------------------
        %qTmp_fac1 = M3_regular( M3_rectangle(qList3, 1, 50,  1) );
        %qTmp_fac2 = M3_regular( featureList12_mrv(ap,:) );
        %qTmp = M3_addZero(qTmp_fac1, qTmp_fac2)
        %----------------------------------
        %q6_fac1 = M3_regular( M3_rectangle(qTmp, 1, 8,  1) );
        %q6_fac2 = M3_regular( featureList2_qd(ap,:) );
        %qList6 = M3_addZero(q6_fac1, q6_fac2);
        %q6Num = sum(qList6 >0)


        %q7_fac1 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, 50,  1) );
        %q7_fac2 = M3_regular( featureTmp2 );
        %qList7 = M3_addZeroZ([q7_fac1; q7_fac2]);
        %----------------------------------
        %qList7 = featureList1_qd(ap  ,:);
        %----------------------------------
        %q7_fac1 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, 50,  1) );
        %q7_fac2 = M3_regular( featureTmp1 );
        %qList7 = M3_addZero(q7_fac1, q7_fac2)
        %----------------------------------
        %q7_fac1 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, 50,  1) );
        %q7_fac2 = M3_regular( featureTmp3 );
        %----------------------------------
        %q7_fac1 = M3_regular( M3_rectangle(featureList0_qd(ap,:), 1, 10,  1) );
        %q7_fac2 = M3_regular( featureList9_mcao(ap,:) )
        %qList7 = M3_addZero(q7_fac1, q7_fac2);
        %q7Num = sum(qList7 >0)
        %----------------------------------
        %q7_fac1 = M3_rectangle(featureList7_cje(ap,:), 1, 100,  1);
        %q7_fac1 = M3_rectangle(featureList0_qd(ap,:), 1, 10, 1);
        %q7_fac2 = 1  - M3_regular( featureList10_activity(ap,:) );
        %qList7 = M3_addZeroZ([q7_fac1; q7_fac2])
        %q7Num = sum(qList7 >0)
        %----------------------------------
        %qList7 = featureList5_hcl(ap, :);
        %----------------------------------
        %qList7_fac1 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, 300,  1) );
        %qList7_fac2 = M3_regular( featureList5_hcl(ap,:) );
        %qList7_fac1 = M3_regular( M3_rectangle(featureTmp1, 1, 30,  1) );
        %q7_fac1 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, 100,  1) );
        %q7_fac1 = M3_regular( M3_rectangle(featureList5_hcl(ap,:), 1, 20,  1) );
        %q7_fac2 = M3_regular( featureList2_qd(ap,:) );
        %q7_fac1 = M3_regular( M3_rectangle(qList5, 1, 30,  1) );
        %q7_fac1 = M3_regular( M3_rectangle(qList5, 1, numHots,  1) );
        %q7_fac1 = M3_regular( M3_rectangle(qList5, 1, 20,  1) );
        %q7_fac2 = 1 - M3_regular( featureList11_cjeRal(ap,:) );
        %q7_fac2 =  M3_regular( featureList11_cjeRal(ap,:) );
        %q7_fac2 = M3_regular( featureList2_qd(ap,:) );
        %qList7 = M3_addZero(q7_fac1, q7_fac2)
        %q7Num = sum(qList7 >0)
        %----------------------------------
        %q7_fac1 = M3_regular( M3_rectangle(qList3, 1, 30,  1) );
        %q7_fac2 = M3_regular( featureList3(ap,:) );
        %qList7 = M3_addZero(q7_fac1, q7_fac2);
        %q7Num = sum(qList7 >0)
        %q7_fac1 = M3_rectangle(featureTmp1,            1, 600, 1);
        %q7_fac2 = M3_rectangle(featureList7_cje(ap,:), 1, 500, 1);
        %q7_fac1 =  M3_rectangle(qList4, 1, 6,  1);
        %q7_fac2 = M3_regular( featureList7_cje(ap,:) );
        %q7_fac1 = M3_regular( M3_rectangle(qTmp, 0, 8,  1) );
        %q7_fac4 = M3_regular( featureList2_qd(ap,:) );
        %q7_fac2 =  M3_regular( featureList5_hcl(ap,:) );
        %q7_fac2 =  1 - M3_regular( featureList3(ap,:) );
        %q7_fac4 = M3_regular(featureList2_qd(ap,:));
        %q7_fac2 =  M3_regular( featureList1_qd(ap,:) );
        %q7_fac2 = M3_regular( featureList11_cjeRal(ap,:) );
        %q7_fac2 =  M3_regular( featureList7_cje(ap,:) );
        %qList7 = M3_addZeroZ([q7_fac1; q7_fac2; q7_fac4]);
        %q7Num = sum(qList7 >0)
        %----------------------------------
        q7_fac1 = M3_rectangle(qList4, 1, 3,  1);
        q7_fac2 = M3_regular( featureList3(ap,:) );
        qList7 = M3_addZero(q7_fac1, q7_fac2);
        q7Num = sum(qList7 >0)





        %if (ratioPositive > 0.4)
        %    qList8 = qList3;
        %else
        %    qList8 = qList7;
        %end
        %----------------------------------
        %qList8 =featureList8_cjeRal(ap,:);
        %iidx = M3_getFilter(featureList7_cje(ap,:), 1, 50,   1);
        %q8Num_filter1 = length(iidx)
        %qList8(iidx) = 0;
        %q8Num = sum(qList8 >0)
        %----------------------------------
        %q8_fac1 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, 50,  1) );
        %q8_fac2 = M3_regular( featureTmp1 );
        %q8_fac3 = M3_regular( featureList1_qd(ap,:) );
        %qList8 = M3_addZeroZ([q8_fac1; q8_fac2; 0.5*(q8_fac3)]);
        %----------------------------------
        %q8_fac1 = M3_regular( M3_rectangle(featureList0_qd(ap,:), 1, 10,  1) );
        %q8_fac2 = M3_regular( featureList7_cje(ap,:) )
        %qList8 = M3_addZero(q8_fac1, q8_fac2);
        %----------------------------------
        %qList8 = featureList2_qd(ap, :);
        %q8Num = sum(qList8 >0)
        %----------------------------------
        %q8_fac1 = M3_rectangle(featureList0_qd(ap,:), 1, 100, 1) ;
        %q8_fac2 = 1  - M3_regular( featureList10_activity(ap,:) );
        %qList8 = M3_addZeroZ([q8_fac1; q8_fac2]);
        %q8Num = sum(qList8 >0)
        %----------------------------------
        %q8_fac1 = M3_regular( M3_rectangle(featureList0_qd(ap,:), 1, 10, 1) );
        %q8_fac2 = M3_regular( featureList9_mcao(ap,:) );
        %qList8 = M3_addZeroZ([q8_fac1; q8_fac2]);
        %q8Num = sum(qList8 >0)
        %----------------------------------
        %q8_fac1 = M3_regular( M3_rectangle(featureList5_hcl(ap,:), 1, 100,  1) );
        %q8_fac1 = M3_regular( M3_rectangle(qList3, 1, 15,  1) );
        %q8_fac1 = M3_regular( M3_rectangle(qList5, 1, 12,  1) );
        %q8_fac2 = 1 - M3_regular( featureList12_mrv(ap,:) );
        %qList8 = M3_addZero(q8_fac1, q8_fac2)
        %q8Num = sum(qList8 >0)
        %----------------------------------
        %qList8_fac0 = M3_regular( M3_rectangle(featureList7_cje(ap,:), 1, 500,  1) );
        %qList8_fac1 = M3_rectangle_mannul(featureList2_qd(ap,:));
        %qList8_fac2 = M3_regular( M3_rectangle(featureTmp1, 1, 1200,  1) );
        %qList8_fac3 = M3_regular( featureList2_qd(ap,:) );
        %qList8 = M3_addZeroZ([ qList8_fac0; qList8_fac1; qList8_fac2; qList8_fac3]);
        %q8Num = sum(qList8 >0)
        %----------------------------------
        %q8_fac1 = M3_regular( M3_rectangle(qList2, 1, 3,  1) );
        %q8_fac2 = M3_regular( featureList6_x(ap,:) );
        %q8_fac2 = M3_regular(q_cje_normal);
        %qList8 = M3_addZeroZ([ q8_fac1; q8_fac2; ]);
        %q8Num = sum(qList8 >0)
        %----------------------------------
        q6_fac1 = M3_rectangle(qList5, 1, 3,  1);
        q6_fac2 = M3_regular( featureList3(ap,:) );
        qList6 = M3_addZero(q6_fac1, q6_fac2);
        q6Num = sum(qList6 >0)

        % step from q8
        %q6_fac1 = M3_regular( M3_rectangle(qList4, 1, 3,  1) );
        %q6_fac1 = M3_regular( M3_rectangle(qList5, 1, 6,  1) );
        %q6_fac1 = M3_regular( M3_rectangle(qList2, 1, 20,  1) );
        %q6_fac2 = M3_regular( featureList2_qd(ap,:) );
        %q6_fac2 = 1 - M3_regular( featureList12_mrv(ap,:) );
        %q6_fac2 =  1 - M3_regular( featureTmp1 );
        %q6_fac1 = M3_regular( M3_rectangle(qList2, 1, 3,  1) );
        %q6_fac2 = 1 - M3_regular( featureList12_mrv(ap,:) );
        %q6_fac2 =  M3_regular( featureList6_x(ap,:) );
        %q6_fac2 = 1 - M3_regular( featureList11_cjeRal(ap,:) );
        %q6_fac2 = M3_regular( featureList7_cje(ap,:) );
        %qList6 = M3_addZero(q6_fac1, q6_fac2);
        %q6Num = sum(qList6 >0)
        %----------------------------------
        %qList6 = qList2;
        %qList6 = featureList6_x(ap,:);
        %----------------------------------
        %q8_fac1 = M3_rectangle(featureList3(ap,:), 1, 100,  1);
        %q8_fac2 = M3_regular( featureList2_qd(ap,:) );
        %qList8  = M3_addZeroZ([ q8_fac1; q8_fac2]);
        %q8Num = sum(qList8 >0)
        %----------------------------------
        %q8_fac1 = M3_rectangle(qList5, 1, 3,  1);
        %q8_fac2 = M3_regular( featureList3(ap,:) );
        %qList8 = M3_addZero(q8_fac1, q8_fac2);
        %q8Num = sum(qList8 >0)
        %q8_fac1 = 1 - M3_rectangle(featureTmp5, 1, 300,  1);
        q8_fac1 = M3_rectangle(featureTmp5, 1, 300, 0);
        q8_fac2 = M3_regular(featureList2_qd(ap,:));
        qList8 = M3_addZeroZ([q8_fac1; q8_fac2]);
        q8Num = sum(qList8 >0)


        q2_fac1 = M3_rectangle(qList8, 1, 3,  1);
        q2_fac2 = M3_regular( featureList3(ap,:) );
        qList2 = M3_addZero(q2_fac1, q2_fac2);
        q2Num = sum(qList2 >0)




        % ******************************************
        %   market indice
        % ******************************************
        qtmp_=sort(qList2, "descend");
        marketIndi1=qtmp_(1:10)
        marketIndi1_s=std(marketIndi1)
        marketIndi1_m=mean(marketIndi1)
        printf("AP %2d, marketindi(m:%f \t s:%f)\n", ap, marketIndi1_m, marketIndi1_s);


        qList1 = rand(size(qList4));
        %qList3 = qList4;
        %qList8 = qList2;

        % ******************************************
        %   record valid count.
        % ******************************************
        checkZoneUnitWidth = 30;
        baseCol = indexEnd + checkZoneUnitWidth*(0) + 23 +7;
        lastCol = indexEnd + checkZoneUnitWidth*(7) + 23 +7;
        validArray = [ sum(qList1>0)   sum(qList2>0)   sum(qList3>0)   sum(qList4>0)   sum(qList5>0)   sum(qList6>0)   sum(qList7>0)   sum(qList8>0) ]
        m3(ap, baseCol:checkZoneUnitWidth:lastCol ) = validArray;
        %m3(ap, baseCol + 5 ) = 123.123;
        %m3(ap, baseCol + 6 ) = sucRate;
        %m3(ap, baseCol + 7 ) = sum(qList1>0);
        %if (ratioPositive < 0.3 || roiRatio < -0.02 || lastRoiRatio < -0.02)
        if (roiRatio < -0.015 || lastRoiRatio < -0.015)
            validArray(2) = 0;
            validArray(3) = 0;
            validArray(4) = 0;
            validArray(5) = 0;
            validArray(6) = 0;
            validArray(7) = 0;
            validArray(8) = 0;
            m3(ap, baseCol:checkZoneUnitWidth:lastCol ) = validArray;
            %m3(ap, baseCol:checkZoneUnitWidth:lastCol) = zeros(size(validArray));
        end
        



        % ******************************************
        %   copy the model3-sort data to somewhere.
        % ******************************************
        if (0)
            % ****************
            %  1.
            % ****************
            needRankList = qList;
            n = length(needRankList);
            [a I] = sort(needRankList);
            tmp(I(:)) = 1:n;
            sortRankList = tmp/m3_ids;
    
            indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue - 0;
            indexL = indexStart:m3_id_width:indexEnd ;
            m3(ap, indexL) =  sortRankList;


            % ****************
            %   2. roi - rank
            % ****************
            indexStart  = m3_header_width + m3_index_strategy1_roi;
            indexL = indexStart:m3_id_width:indexEnd ;
            roiList = m3(ap, indexL);
       
            needRankList = roiList;
            n = length(needRankList);
            [a I] = sort(needRankList);
            tmp(I(:)) = 1:n;
            sortRankList = tmp/m3_ids;


            indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue - 2;
            indexL = indexStart:m3_id_width:indexEnd ;
            m3(ap, indexL) =  sortRankList;
        end


        % ******************************************
        %  save to y-area.(not sort)
        % ******************************************
        iModel = 2;
        indexStart =  m3_header_width + m3_index_strategy2_b0 + iModel ;
        indexL = indexStart:m3_id_width:indexEnd ;
        m3(ap, indexL) = qList2;

        iModel = 3;
        indexStart =  m3_header_width + m3_index_strategy2_b0 + iModel ;
        indexL = indexStart:m3_id_width:indexEnd ;
        m3(ap, indexL) = qList3;

        iModel = 4;
        indexStart =  m3_header_width + m3_index_strategy2_b0 + iModel ;
        indexL = indexStart:m3_id_width:indexEnd ;
        m3(ap, indexL) = qList4;

        iModel = 5;
        indexStart =  m3_header_width + m3_index_strategy2_b0 + iModel ;
        indexL = indexStart:m3_id_width:indexEnd ;
        m3(ap, indexL) = qList5;

        iModel = 6;
        indexStart =  m3_header_width + m3_index_strategy2_b0 + iModel ;
        indexL = indexStart:m3_id_width:indexEnd ;
        m3(ap, indexL) = qList6;

        iModel = 7;
        indexStart =  m3_header_width + m3_index_strategy2_b0 + iModel ;
        indexL = indexStart:m3_id_width:indexEnd ;
        m3(ap, indexL) = qList7;

        iModel = 8;
        indexStart =  m3_header_width + m3_index_strategy2_b0 + iModel ;
        indexL = indexStart:m3_id_width:indexEnd ;
        m3(ap, indexL) = qList8;



        % ****************************************************
        %  update valid column.
        % ****************************************************
        m3(ap, m3_index_global_valid)  =  countValid;


        % ****************************************************
        %  at once, sort this.
        %   reason:
        %    1. avoid big number of qList in the m3.edat.
        % ****************************************************
        iModel = 2;
        rankList = qList2;
        n = length(rankList);
        [a I] = sort(rankList);
        printf("qList2 best3: %d,\t %d,\t %d\n", I(n), I(n-1), I(n-2));
        tmp(I(:)) = 1:n;
        sortRankList = tmp ./ n;
        %indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue + iModel ;
        indexStart =  m3_header_width + m3_index_strategy2_b0+ iModel ;
        indexL = indexStart:m3_id_width:indexEnd ;
        m3(ap, indexL) = sortRankList;


        iModel = 3;
        rankList = qList3;
        n = length(rankList);
        [a I] = sort(rankList);
        tmp(I(:)) = 1:n;
        sortRankList = tmp ./ n;
        %indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue + iModel ;
        indexStart =  m3_header_width + m3_index_strategy2_b0+ iModel ;
        indexL = indexStart:m3_id_width:indexEnd ;
        m3(ap, indexL) = sortRankList;

        iModel = 4;
        rankList = qList4;
        n = length(rankList);
        [a I] = sort(rankList);
        tmp(I(:)) = 1:n;
        sortRankList = tmp ./ n;
        %indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue + iModel ;
        indexStart =  m3_header_width + m3_index_strategy2_b0+ iModel ;
        indexL = indexStart:m3_id_width:indexEnd ;
        m3(ap, indexL) = sortRankList;

        iModel = 5;
        rankList = qList5;
        n = length(rankList);
        [a I] = sort(rankList);
        tmp(I(:)) = 1:n;
        sortRankList = tmp ./ n;
        %indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue + iModel ;
        indexStart =  m3_header_width + m3_index_strategy2_b0+ iModel ;
        indexL = indexStart:m3_id_width:indexEnd ;
        m3(ap, indexL) = sortRankList;

        iModel = 6;
        rankList = qList6;
        n = length(rankList);
        [a I] = sort(rankList);
        tmp(I(:)) = 1:n;
        sortRankList = tmp ./ n;
        %indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue + iModel ;
        indexStart =  m3_header_width + m3_index_strategy2_b0+ iModel ;
        indexL = indexStart:m3_id_width:indexEnd ;
        m3(ap, indexL) = sortRankList;

        iModel = 7;
        rankList = qList7;
        n = length(rankList);
        [a I] = sort(rankList);
        printf("qList7 best3: %d,\t %d,\t %d\n", I(n), I(n-1), I(n-2));
        tmp(I(:)) = 1:n;
        sortRankList = tmp ./ n;
        %indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue + iModel ;
        indexStart =  m3_header_width + m3_index_strategy2_b0+ iModel ;
        indexL = indexStart:m3_id_width:indexEnd ;
        m3(ap, indexL) = sortRankList;

        iModel = 8;
        rankList = qList8;
        n = length(rankList);
        [a I] = sort(rankList);
        printf("qList8 best3: %d,\t %d,\t %d\n", I(n), I(n-1), I(n-2));
        tmp(I(:)) = 1:n;
        sortRankList = tmp ./ n;
        %indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue + iModel ;
        indexStart =  m3_header_width + m3_index_strategy2_b0+ iModel ;
        indexL = indexStart:m3_id_width:indexEnd ;
        m3(ap, indexL) = sortRankList;



    end  % endof ap


    end  % endof beta.


    % ********************
    %  Note
    % ********************
end

endfunction


######################################################## FILE.
function  [] = M3_makeR1R2R3()
######################################################## FILE.
global m3 ;
global m3_ids ;
global nGroups;
global m2_groupSize;

global m3_header_width;
global m3_index_strategy2_y_estimate;
global m3_index_strategy2_y_estimate2;
global m3_index_strategy1_roi;
global m3_id_width;
global m3_index_date1     ;    % header area
global m3_index_date2     ;    % header area
global m3_index_global_valid;  % header area
global m3_index_envStress ;%= 3;    % header area
global m3_index_r1        ;%= 5;   % header area

global m3_index_strategy2_b1;
global m3_index_strategy2_b0;
global m3_index_strategy2_b0_pvalue;
global m3_index_strategy2_y ;
global m3_index_strategy1_y ;

global m2;
global m2_header_width;
global m2_index_p;
global m2_id_width;

global stockIdStr;
global G_pg;
global G_ModelSelect_toFile;
global G_MaxModels ;
global A_table;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% make data structure 'yyy' for statistic info.
%% 
%% 
%%  Predict-Unit:
%%  --------------------
%%   1. (ap, sid)
%%   2. with 6 model predict results.
%%   3. with real-y, that is not in terms of roi.
%%   4. (real-y, y-roi)  co-inherence relation!
%% 
%% 
%%  --------------------------------------------------
%% 
%% 0    : y_est
%% 1    : y_model (sync)
%% 2    : y_real
%% 3-8  : y model-x (already regurlize)
%% 9-10 : reserve.
%% 11   : y -roi. ( co-inherence with y_real )
%% 12   : ap
%% 13   : sid
%% 14   : roi-future.
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyItems=0;
yyy = zeros(1,15);
for ap=1:nGroups    % ap means!


    % *******************
    %  validation.
    % *******************
    if (m3(ap, m3_index_global_valid) == 0)
        continue;
    end


    % *******************
    %  fill data to yyy[].
    % *******************
    for sid=1:m3_ids
       indexStart =  m3_header_width + m3_index_strategy2_y_estimate + (sid-1)*m3_id_width;
       item = m3(ap, indexStart+(0:10));

       % ******************************************************
       %  model [1,6] come from m3_index_strategy2_b0_pvalue
       % ******************************************************
       indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue + (sid-1)*m3_id_width;
       item(4:11)  =  m3(ap, indexStart+(1:8));



       % ******************************************************
       %   generate roi, and store to m3_index_strategy1_roi
       % ******************************************************
       if (1)
           indexRoi = m3_header_width + m3_index_strategy1_roi + (sid-1)*m3_id_width;
           thisRoi = m3(ap,indexRoi);
           indexY_real = m3_header_width + m3_index_strategy1_y + (sid-1)*m3_id_width;
           Y_real = m3(ap, indexY_real);
           item = [item  thisRoi ap sid 0 ] % 15 elements.
           item(3)   =  Y_real;
           item(15)  =  Y_real;
       end

       % ******************************************************
       %   generate roi, and store to m3_index_strategy1_roi
       %    delete...
       %    delete...
       %    delete...
       %    delete...
       % ******************************************************


       % means est>0, model>0, real>0.
       %if (item(1)>0 && item(2)>0 && item(3)>0)
       if (1>0)
           yyItems++;
           item
           yyy(yyItems, :) = item;
       end
    end
end

if (0)
numUnits = size(yyy, 1);
for u=1:numUnits
    ap  = yyy(u, [13]);
    sid = yyy(u, [14]);

    if (ap+1 > nGroups)
        continue;
    end
    % *******************
    %  fill 
    % *******************
    sid
    indexRoi = m3_header_width + m3_index_strategy1_roi + (sid-1)*m3_id_width
    yyy(u, 15) = m3(ap+1,indexRoi);  % update roi-future.

end
end



%yyy(:,1) = yyy(:,1)/m3_ids;  % y est
%yyy(:,3) = yyy(:,3)/m3_ids;  % y real
format bank;
sizeYYY=size(yyy)
format short;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Pg calculate.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rowNum=1;
%for const_rank=[0.8 0.9 0.95]
for const_rank=[0.990 0.995 0.999]
    
    idxx = yyy(:, 4) > const_rank;  % model 1
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY>0.8) ;
    s2 = sum(idxx);
    Pg(rowNum, 1) = s1 / s2;
    
    idxx = yyy(:, 5) > const_rank;  % model 2
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY>0.8);
    s2 = sum(idxx);
    Pg(rowNum, 2) = s1 / s2;
     
    idxx = yyy(:, 6) > const_rank;  % model 3
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY>0.8) ;
    s2 = sum(idxx);
    Pg(rowNum, 3) = s1 / s2;
       
     
    idxx = yyy(:, 7) > const_rank;  % model 4
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY>0.8) ;
    s2 = sum(idxx);
    Pg(rowNum, 4) = s1 / s2;
    
     
    idxx = yyy(:, 8) > const_rank;  % model 5
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY>0.8) ;
    s2 = sum(idxx);
    Pg(rowNum, 5) = s1 / s2;
    
    
    idxx = yyy(:, 9) > const_rank;  % model 6
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY>0.8) ;
    s2 = sum(idxx);
    Pg(rowNum, 6) = s1 / s2;


    idxx = yyy(:, 10) > const_rank;  % model 7
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY>0.8) ;
    s2 = sum(idxx);
    Pg(rowNum, 7) = s1 / s2;


    idxx = yyy(:, 11) > const_rank;  % model 8
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY>0.8) ;
    s2 = sum(idxx);
    Pg(rowNum, 8) = s1 / s2;

    
    rowNum++;
end


testCounts = size(yyy)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Pb calculate.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rowNum=1;
%for const_rank=[0.8 0.9 0.95]
for const_rank=[0.990 0.995 0.999]
    
    idxx = yyy(:, 4) > const_rank;  % model 1
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY<0.5) ;
    s2 = sum(idxx);
    Pb(rowNum, 1) = s1 / s2;
    
    idxx = yyy(:, 5) > const_rank;  % model 2
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY<0.5) ;
    s2 = sum(idxx);
    Pb(rowNum, 2) = s1 / s2;
     
    idxx = yyy(:, 6) > const_rank;  % model 3
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY<0.5) ;
    s2 = sum(idxx);
    Pb(rowNum, 3) = s1 / s2;
       
     
    idxx = yyy(:, 7) > const_rank;  % model 4
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY<0.5) ;
    s2 = sum(idxx);
    Pb(rowNum, 4) = s1 / s2;
    
     
    idxx = yyy(:, 8) > const_rank;  % model 5
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY<0.5) ;
    s2 = sum(idxx);
    Pb(rowNum, 5) = s1 / s2;
    
    
    idxx = yyy(:, 9) > const_rank;  % model 6
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY<0.5) ;
    s2 = sum(idxx);
    Pb(rowNum, 6) = s1 / s2;


    idxx = yyy(:, 10) > const_rank;  % model 7
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY<0.5) ;
    s2 = sum(idxx);
    Pb(rowNum, 7) = s1 / s2;


    idxx = yyy(:, 11) > const_rank;  % model 8
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY<0.5) ;
    s2 = sum(idxx);
    Pb(rowNum, 8) = s1 / s2;
    
    
    rowNum++;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Pg_ext calculate.
%%
%%   the real model point 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rowNum=1;
%for const_rank=[0.8 0.9 0.95 0.98 0.99]
for const_rank=[0.990 0.995 0.999]
    
    idxx = yyy(:, 4) > const_rank;  % model 1
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY) ;
    s2 = sum(idxx)
    Pg_ext(rowNum, 1) = s1 / s2;
    
    idxx = yyy(:, 5) > const_rank;  % model 2
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY);
    s2 = sum(idxx);
    Pg_ext(rowNum, 2) = s1 / s2;
     
    idxx = yyy(:, 6) > const_rank;  % model 3
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY) ;
    s2 = sum(idxx);
    Pg_ext(rowNum, 3) = s1 / s2;
       
    idxx = yyy(:, 7) > const_rank;  % model 4
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY) ;
    s2 = sum(idxx);
    Pg_ext(rowNum, 4) = s1 / s2;
     
    idxx = yyy(:, 8) > const_rank;  % model 5
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY) ;
    s2 = sum(idxx);
    Pg_ext(rowNum, 5) = s1 / s2;
    
    idxx = yyy(:, 9) > const_rank;  % model 6
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY) ;
    s2 = sum(idxx);
    Pg_ext(rowNum, 6) = s1 / s2;

 
    idxx = yyy(:, 10) > const_rank;  % model 7
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY) ;
    s2 = sum(idxx);
    Pg_ext(rowNum, 7) = s1 / s2;


    idxx = yyy(:, 11) > const_rank;  % model 8
    matchY   =  yyy(idxx, 3);
    s1 = sum(matchY) ;
    s2 = sum(idxx);
    Pg_ext(rowNum, 8) = s1 / s2;


    rowNum++;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ROE --> ROI map.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rowNum=1;
%interval = 0.01;
interval = 0.001;
for roe=[0.99:interval:1]

roeHH = roe

    predict = yyy(:,4);  % model 1
    idxx = (predict<=roe & predict>=(roe-interval) );
    matchY   =  yyy(idxx, 15);
    s1 = sum(matchY);
    s2 = sum(idxx)
    e2i(rowNum, 1) = s1 / s2;
    
    predict = yyy(:,5);  % model 2
    idxx = (predict<=roe & predict>=(roe-interval) );
    matchY   =  yyy(idxx, 15);
    s1 = sum(matchY);
    s2 = sum(idxx)
    e2i(rowNum, 2) = s1 / s2;
 
    predict = yyy(:,6);  % model 3
    idxx = (predict<=roe & predict>=(roe-interval) );
    matchY   =  yyy(idxx, 15);
    %info   =  yyy(idxx, [13 14 6 12 15])
    s1 = sum(matchY);
    s2 = sum(idxx)
    e2i(rowNum, 3) = s1 / s2;
 
    predict = yyy(:,7);  % model 4
    idxx = (predict<=roe & predict>=(roe-interval) );
    matchY   =  yyy(idxx, 15);
    s1 = sum(matchY);
    s2 = sum(idxx)
    e2i(rowNum, 4) = s1 / s2;
 
    predict = yyy(:,8);  % model 5
    idxx = (predict<=roe & predict>=(roe-interval) );
    matchY   =  yyy(idxx, 15);
    s1 = sum(matchY);
    s2 = sum(idxx)
    e2i(rowNum, 5) = s1 / s2;
 
 
    predict = yyy(:,9);  % model 6
    idxx = (predict<=roe & predict>=(roe-interval) );
    matchY   =  yyy(idxx, 15);
    s1 = sum(matchY);
    s2 = sum(idxx)
    e2i(rowNum, 6) = s1 / s2;

    predict = yyy(:,10);  % model 7
    idxx = (predict<=roe & predict>=(roe-interval) );
    matchY   =  yyy(idxx, 15);
    s1 = sum(matchY);
    s2 = sum(idxx)
    e2i(rowNum, 7) = s1 / s2;
 
    predict = yyy(:,11);  % model 8
    idxx = (predict<=roe & predict>=(roe-interval) );
    matchY   =  yyy(idxx, 15);
    s1 = sum(matchY);
    s2 = sum(idxx)
    e2i(rowNum, 8) = s1 / s2;
   
fff=333333333
    rowNum++;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ROE --> ROI map.
%%
%%  Objective:
%%    Judge the model whether following or close to the real ROI.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rowNum=1;
for roe=[0.1:0.1:1]

    % *******************
    %  real-y  ---> ROI
    %  Note:
    %    be carefull distinct from stra1-real-y and stra2-real-y.
    % *******************
    predict = yyy(:,3); % real y
    % predict = yyy(:,13); % real y
    idxx = (predict<=roe & predict>=(roe-0.1) );
    matchY   =  yyy(idxx, 15);
    s1 = sum(matchY)
    s2 = sum(idxx)

    e2i_global(rowNum, 1) = s1 / s2;
    e2i_global(rowNum, 8) = 0.88;

    roiList(rowNum) = s1/s2;
    roiNums(rowNum) = s2;
    rowNum++;
end


rowNum=1;
for roe=[0.1:0.1:1]

    % *******************
    %  roi     ---> ROI
    % *******************
    predict = yyy(:,12); % real y
    [a idx] = sort(predict);

    begin = floor((roe-0.1)*length(predict));    
    enddd = floor((roe    )*length(predict));    
    if (begin<1) begin = 1; end;

    e2i_global(rowNum, 2) = sum(a(begin:enddd)) / (enddd-begin+1);
    e2i_global(rowNum, 8) = 0.88;


    roiList2(rowNum) = e2i_global(rowNum, 2);
    roiNums2(rowNum) = enddd-begin+1;
    rowNum++;

end


meanRoi  = mean(yyy(:,12))
meanRoi1 = mean(roiList)
meanRoi2 = mean(roiList2)
roiNums
roiNums2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Proi calculate.
% 
%   1. Do NOT modify Y-est!
%   2. Be able to specify which model to predict!
%
%
%
%  check Zone structure:
% ------------------------------------------------------- 
%  | 9 best ones (predict) | 2 verify ones | sum data |
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Proi = zeros(1,8);
valid_startAp = 0;
validAps=0;
checkZoneUnitWidth=30;

indexEnd   =  m3_header_width + m3_ids*m3_id_width;
%m3(:, indexEnd + (20:59) ) = zeros(nGroups, 40);  % clear.
if (1)
for Modelx = 1:G_MaxModels

sumRoi_aps=0;
sumRoe_aps=0;
sumSuc_aps=0;
validAps=0;
tradeAps=0;
aaa=0;
bbb=0;
for ap=1:nGroups    % ap means!
%for ap=1:(nGroups-1)    % ap means!
            
            % *******************
            %  validation.
            % *******************
            if (m3(ap, m3_index_global_valid) == 0)
                continue;
            end
            if (ap < 2) % at least, 3
                continue;
            end


            % *******************
            %  validation pass.
            % *******************
            validAps ++;
            if ( valid_startAp == 0)
                 valid_startAp = ap;
            end


            % ****************************
            %   support dynamic model!
            % ****************************
            %if (Modelx == 8)
            if (0)
                lastAp = ap - 1;
                lastAp = ap ;  %%%%%%%%%%%%%%%%%%%%% --------------------------- NOTE!

                candidate1 = 2;
                baseCol1 = indexEnd + checkZoneUnitWidth*(candidate1-1) + 23 + 1;
                cand1Indi1 = m3(ap  , baseCol1 ) ;
                cand1Indi2 = m3(ap-1, baseCol1 ) ;

                candidate2 = 7; %3;
                baseCol2 = indexEnd + checkZoneUnitWidth*(candidate2-1) + 23 + 1;
                cand2Indi1 = m3(ap  , baseCol2 ) ;
                cand2Indi2 = m3(ap-1, baseCol2 ) ;


                % *****************
                %   dynamic here.
                % *****************
                if (cand1Indi1 > cand2Indi1 && cand1Indi2>cand2Indi2)
                    chooseModelx = candidate1;
                else
                    chooseModelx = candidate2;
                end

                ap
                chooseModelx


                % *******************
                %  src -> dest 
                % *******************
                indexStart    =  m3_header_width + m3_index_strategy2_b0_pvalue + chooseModelx;
                YEstIdx_y     =  indexStart : m3_id_width : indexEnd; 
                YEstVector_y  =  m3(ap, YEstIdx_y);

                indexStart     =  m3_header_width + m3_index_strategy2_b0_pvalue + 8;
                YEstIdx_yDest  =  indexStart : m3_id_width : indexEnd;
                m3(ap, YEstIdx_yDest) = YEstVector_y;

                % *******************
                %  recover.
                % *******************
            end


            % *************************************
            %   special for dynamic model.
            %     write to M3 for use in future.
            % *************************************
            if (Modelx == 8)
                %indexStart_src =  m3_header_width + m3_index_strategy2_b0_pvalue + chooseModelx;
                %YEstIdx_src    =  indexStart_src : m3_id_width : indexEnd; 
                %indexStart_dest=  m3_header_width + m3_index_strategy2_b0_pvalue + 8;
                %YEstIdx_dest   =  indexStart_dest : m3_id_width : indexEnd; 
                %m3(ap, YEstIdx_dest) = m3(ap, YEstIdx_src);
            end
 


            % ************
            %  make y1
            % ************
            indexStart     =  m3_header_width + m3_index_strategy2_b0_pvalue + Modelx;
            YEstIdx_y1     =  indexStart : m3_id_width : indexEnd; 
            YEstVector_y1  =  m3(ap, YEstIdx_y1);
            


            % *******************
            % update Y-est!
            % *******************
            [a ix] = sort(YEstVector_y1);
            tmp(ix) = 1:m3_ids;
            %m3(ap, YEstIdx) = tmp;


            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% append r1-r9
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            m3(ap, indexEnd + (Modelx-1)*checkZoneUnitWidth + 1) = Modelx + 0.271828;  % start flag of Unit-model.
            bestNums = 9;  %  m3_ids*0.05;
            for k=1:bestNums
                i = ix(m3_ids - k + 1);
                m3(ap, indexEnd + (Modelx-1)*checkZoneUnitWidth + 1 + k) = k*10000 + i + stockIdStr(i) / 1000000;
            end



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% make roi values for r.95 (means r greater than 0.95) 
            %% put in the test-area (testZone).
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            bestNums = floor(m3_ids*0.05);
            bestNums = floor(m3_ids*0.02);
            bestNums = min(2, bestNums); %  maybe useful.
            bestNums = 2; % fix number. 
            bestNums = 3; % fix number. 
            lastAp = ap - 1;



            % ****************************
            % calculate yy for each sid.
            % ****************************
            YY = 0;
            for sid=1:m3_ids

                indexRoi = m3_header_width + m3_index_strategy1_roi + (sid-1)*m3_id_width;
                thisRoi  = m3(ap,indexRoi) ;

                YY(sid) = thisRoi;
            end


            % *********************************
            %  make rank of YY. (YY --- container in terms of 'Rank of Roi')
            % *********************************
            [a ix] = sort(YY);
            YYRank(ix) = 1:m3_ids;
            lastAp
            a
            YY_ = sort(YY)
            
            aaa(lastAp,1:8) = a(m3_ids-7:m3_ids);
            bbb(lastAp,1:8) = a(1:8);


            % *********************************
            %  make roe, roi.
            % *********************************
            rankSum = 0;
            roiSum = 0;
            tradeSum = 0;
            predictSum = 0;
            successSum = 0;
            for k=1:bestNums

                baseCol = indexEnd + checkZoneUnitWidth*(Modelx-1) + 1;
                %sid = m3(lastAp, baseCol + k);
                sid = m3(ap, baseCol + k);
                sid = floor(sid - k*10000);
                %sid = m3(lastAp, baseCol + k +2) ;  % mark !!!!!!!!!!!!!
                %sid = floor(sid - (k+2)*10000);      % mark !!!!!!!!!!!!!

                if (sid<1)
                    continue;
                end


                % if not valid item, do not make ROI.
                validNum = m3(ap, baseCol + 29);
                baseCol = indexEnd + checkZoneUnitWidth*(Modelx-1) + 11 + (k-1)*4;
                if (k>validNum)
                    printf("mask %d\t %d\t roi:%f\n", ap, stockIdStr(sid), YY(sid));

                    m3(ap, baseCol+0) = k ;
                    m3(ap, baseCol+1) = 666777/1000000 ;
                    m3(ap, baseCol+2) = 0.; % thisRoi
                    m3(ap, baseCol+3) = YY(sid);
                
                    %rankSum = rankSum + m3_ids/2;
                    %roiSum  = roiSum + 0;
                    predictSum ++;
                    if (YY(sid) < 0)
                       successSum ++;
                    end
                    continue;
                end


                % stream to file M3.
                m3(ap, baseCol+0) = k ;
                m3(ap, baseCol+1) = stockIdStr(sid)/ 1000000 ;
                m3(ap, baseCol+2) = YY(sid) ; % thisRoi
                %m3(ap, baseCol+3) = YYRank(sid) ;  
                global m3_index_strategy2_indicator1;
                index_time = m3_header_width + m3_index_strategy2_indicator1 + (sid-1)*m3_id_width;
                m3(ap, baseCol+3) = m3(ap,index_time);  

                rankSum = rankSum + YYRank(sid)
                roiSum  = roiSum + YY(sid);

                predictSum ++;
                tradeSum ++;       % do trade, so increase trade number.
                if (YY(sid) > 0)
                   successSum ++;
                end
            end
            %predictSum = sum(m3_ids+1 - (1:bestNums))
            if (tradeSum>0)
                tradeAps ++;
                ROE = rankSum /(tradeSum*m3_ids)  % rank of earn
                ROI = roiSum/tradeSum             % return of invest
            else
                ROE = 0  % rank of earn
                ROI = 0  % return of invest
            end
            sucRate = successSum/predictSum;


            % **********************************************
            %  for each model:
            %    record roi, roe to testZone-header.
            %    record roi-rank
            % **********************************************
            tt = m3(ap, m3_index_date2) ; %
            st = localtime(tt);
            baseCol = indexEnd + checkZoneUnitWidth*(Modelx-1) + 23;
            m3(ap, baseCol + 0 ) = st.wday + 0.001234567 + 0.1*Modelx;
            m3(ap, baseCol + 1 ) = ROE;  % RANK OF EARN.(In other word, roi-rank)
            m3(ap, baseCol + 3 ) = ROI;
            m3(ap, baseCol + 4 ) = sum(YY(:))/length(YY);  % roi-index.
            m3(ap, baseCol + 5 ) = 123.123;
            m3(ap, baseCol + 6 ) = sucRate;
             
             
            % **********************************************
            %   mean roe 
            % **********************************************
            roeListInWindow = m3(ap-2:ap, baseCol + 1 ) ;
            nL = sum(roeListInWindow>0);
            if (nL>0) 
                m3(ap, baseCol + 2 ) = sum(roeListInWindow)/nL;
            else
                m3(ap, baseCol + 2 ) = 0;
            end


            % **********************************************
            %   mean roi (market mean) 
            % **********************************************
            roiListInWindow = m3(ap-2:ap, baseCol + 4 ) ;
            nL = sum(roiListInWindow>0);
            if (nL>0)
                %m3(ap, baseCol + 5 ) = sum(roiListInWindow)/nL;
            else
                %m3(ap, baseCol + 5 ) = 0;
            end



            sumRoe_aps = sumRoe_aps + ROE;
            sumRoi_aps = sumRoi_aps + ROI;
            sumSuc_aps = sumSuc_aps + sucRate;
            A_table(ap,Modelx+0) = ROE;
            %A_table(ap,1)        = ap;

end  % endof each ap.

            % **********************************************
            %  record (mean-) roi, roe to logFile (NOT testZone-header).
            % **********************************************
            col = Modelx;
            %ROEMEAN = sum(m3(:,indexEnd+21))/nGroups
            %ROIMEAN = sum(m3(:,indexEnd+22))/nGroups
            printf("APS: (valid:%d VS trade:%d)\n", validAps, tradeAps);
            if (tradeAps > 0)
                ROEMEAN =  sumRoe_aps/tradeAps;
                ROIMEAN =  sumRoi_aps/tradeAps;
            else
                ROEMEAN =  0;
                ROIMEAN =  0;
            end
            SUCMEAN =  sumSuc_aps/validAps;
            roeCol = indexEnd + checkZoneUnitWidth*(Modelx-1) + 23 + 1;
            Proi(1, col) = ROEMEAN ;
            Proi(2, col) = std(m3(:, roeCol)); % ROE std.
            Proi(3, col) = ROIMEAN ;
            Proi(4, col) = sumRoi_aps;
            Proi(5, col) = SUCMEAN;
            Proi(6, col) = sum(m3(:, indexEnd + checkZoneUnitWidth*(Modelx-1) + 13)); % ROI of champion.

            aaa
            zhuiZhangYiYuan=[ mean(aaa,2) mean(bbb,2) ]

end  % endof each model.
end



% **************************************
%   ....
%   ....
%   ....
% **************************************
maxcols = G_MaxModels;
G_pg = [Pg; zeros(2,maxcols);       ...
        Pg_ext; zeros(2,maxcols);   ...
        Pb; zeros(2,maxcols);       ...
        e2i; zeros(2,maxcols);      ...
        e2i_global; zeros(2,maxcols);  ...
        Proi ];
meanROE = mean(A_table);
mediROE = median(A_table);
A_table = [A_table ; 
           zeros(2, size(A_table,2)); meanROE; mediROE;
];


endfunction




######################################################## FILE.
function  [] = M3_makeChooseSid(policy)
######################################################## FILE.
global m3 ;
global m3_ids ;
global nGroups;
global m2_groupSize;

global m3_header_width;
global m3_index_strategy2_y_estimate;
global m3_id_width;

global m3_index_date1     ;    % header area
global m3_index_date2     ;    % header area


global m3_index_envStress ;%= 3;    % header area
global m3_index_r1        ;%= 5;   % header area

global m3_index_strategy2_b0  ;
global m3_index_strategy2_b1  ;
global m3_index_strategy2_b2  ;
global m3_index_strategy2_b3  ;

global m3_index_global_valid;
global m3_index_strategy2_b0_pvalue  ;

global m2;
global m2_header_width;
global m2_index_p;
global m2_id_width;
global G_lastSid;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear area.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m3(1:nGroups, m3_index_r1) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%  lastAp:    Y_sid... Yest... 
%% 
%%  ap:        Praw
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for specifyRank=1:1                % super parameter
    for iLow =1:1                      % super parameter, always be same as std of neighbor.
    for iUp=1:1                        % super parameter


        for ap=1:nGroups    % ap means!

            %%%%%%%%%%%%%%%%%%%%
            %% m2 timeseries.
            %%%%%%%%%%%%%%%%%%%%
            thisAPidx  = (m2(:,2) == ap);
            timeSeries = m2(thisAPidx, 1); % set timestamp of each ap.
            tt = timeSeries(1);
            st = localtime(tt);
            ymd1 = (st.year+1900)*10000 + (st.mon+1)*100 + st.mday;
            m3(ap, m3_index_date1) = ymd1;
            tt = timeSeries(length(timeSeries));
            st = localtime(tt);
            ymd2 = (st.year+1900)*10000 + (st.mon+1)*100 + st.mday;
            %m3(ap, m3_index_date2) = ymd2 + ceil(st.wday)/10;
            m3(ap, m3_index_date2) = tt; %
            m3(ap, m3_index_date2) = ymd2;%tt; %


            %%%%%%%%%%%%%%%%%%%%
            %% get stock bucket 2: get the last(choose) ones.
            %%%%%%%%%%%%%%%%%%%%
            indexStart =  m3_header_width + m3_index_strategy2_y_estimate; % est
            indexEnd   =  m3_header_width + m3_ids*m3_id_width;
            YEstIdx    =  indexStart : m3_id_width : indexEnd; 
            YEstVector   =  m3(ap, YEstIdx);
            YRealVector  =  m3(ap, YEstIdx-1);
            YModelVector =  m3(ap, YEstIdx+1);


            %%%%%%%%%%%%%%%%%%%%
            %% Validation:
            %% 
            %%  if not pass, then do nothing, and let 'm3_index_r1' be zero!
            %%%%%%%%%%%%%%%%%%%%
            if (m3(ap, m3_index_global_valid)  == 0)
                continue;
            end


            % *********************************************
            % use Y-est-mode-x substitude Y-est!
            %        dynamic select model.
            %       
            %  override YEstVector.
            % *********************************************
            global G_ModelSelect;
            indexStart =  m3_header_width + m3_index_strategy2_b0_pvalue + G_ModelSelect; % est y  (model x)
            indexEnd   =  m3_header_width + m3_ids*m3_id_width;
            YEstIdx_y1     =  indexStart : m3_id_width : indexEnd; 
            YEstVector_y1  =  m3(ap, YEstIdx_y1);

            [a ix] = sort(YEstVector_y1);
            tmp(ix) = 1:m3_ids;
            YEstVector_y1  =  tmp;
            YEstVector     =  YEstVector_y1;  % Note!
            


            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% important!
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            chooseY    =  YEstVector >= (m3_ids*0.95);


            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% record select stock!
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            global stockIdStr;
            stockSelect(ap, :) = (stockIdStr(chooseY))';


            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% cxw policy random:
            %%    select only randon one from good pool.
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (policy == 1)
                chooseY    =  YEstVector >= (m3_ids*0.6);
                %stockSelect(ap, :) = (stockIdStr(chooseY))';

                nn = length(chooseY);
                numGoods = sum(chooseY)
                randii = 1 + floor(rand(1) * numGoods);
                randii_idx = 0;
                for _i=1:nn
                    if (chooseY(_i) == 1)
                        randii_idx = randii_idx + 1;
                    end
    
                    if (randii_idx == randii)
                        break;
                    end
                end
                chooseY=repmat([0], nn, 1);
                chooseY(_i) = 1;

%                chooseY
%                chooseY = YEstVector == (m3_ids)
%                chooseY = YEstVector > (m3_ids - 3)
%                stockIdStr(chooseY) 
%                %% special delete 600104
%                if (600104 ==  stockIdStr(_i) )
%                
%                end

            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% cxw policy random:
            %%    del randon 2 from good pool.
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (policy == 2)
                nn = length(chooseY);
                for _delOne=1:4
                    numGoods = sum(chooseY)
                    randii = 1 + floor(rand(1) * numGoods);
                    randii_idx = 0;
                    for _i=1:nn
                        if (chooseY(_i) == 1)
                            randii_idx = randii_idx + 1;
                        end
        
                        if (randii_idx == randii)
                            break;
                        end
                    end
                    %chooseY=repmat([0], nn, 1);
                    chooseY(_i) = 0;
                end
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% cxw policy random:
            %%    select more better.
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (policy == 3)
                randii = rand(1) ;

                if (randii > 0.5)
                    chooseY    =  YEstVector == (m3_ids - 0);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%% note
                else
                    chooseY    =  YEstVector == (m3_ids - 1);
                    %chooseY    =  YEstVector == (m3_ids - 5);
                end
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% cxw policy (special for regular):
            %%    
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (policy == 4)
                chooseY    =  YEstVector == (m3_ids);
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% cxw policy (special for regular):
            %%    
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (policy == 5)
                %chooseSid1    =  YEstVector == (m3_ids);
                %chooseSid2    =  YEstVector == (m3_ids-1);

                randii = rand(1) ;
                if (randii > 0.5)
                    marginVar = 1;
                else
                    marginVar = 2;
                end

                [_a idx] = sort( YEstVector == (m3_ids) );
                chooseSid1 = idx(m3_ids);
                [_a idx] = sort( YEstVector == (m3_ids-marginVar) );
                chooseSid2 = idx(m3_ids);
               
                PrePreAp = ap - 2 ;   %% important!!!
                indexYest1 =  m3_header_width + m3_id_width*(chooseSid1-1) + m3_index_strategy2_y_estimate;
                indexYest2 =  m3_header_width + m3_id_width*(chooseSid2-1) + m3_index_strategy2_y_estimate;
                indexYRel1 =  indexYest1 - 1;
                indexYRel2 =  indexYest2 - 1;
                   
                delta1 = abs(m3(PrePreAp, indexYRel1) - m3_ids); 
                delta2 = abs(m3(PrePreAp, indexYRel2) - m3_ids + 1);

                if (delta1 < delta2 )
                    chooseY = YEstVector == (m3_ids);
                else
                    chooseY = YEstVector == (m3_ids-marginVar);
                end
            end



            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% cxw policy ( simulate real world):
            %%    
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (policy == 6)
                chooseY    =  YEstVector >= (m3_ids)*0.95;
                %chooseY    =  YEstVector >= (m3_ids)*0.90;

                %chooseY1    =  YEstVector == (m3_ids - 0);
                %chooseY2    =  YEstVector == (m3_ids - 1);
                %chooseY3    =  YEstVector == (m3_ids - 2);
                %chooseY = (chooseY1 | chooseY2 | chooseY3);


                foundHard = 0;
                for sid=1:m3_ids
                    % **************
                    % in the pool.
                    % **************
                    if (chooseY(sid) && (G_lastSid == sid))
                        foundHard = 1;
                        chooseY = zeros(size(chooseY));
                        chooseY(sid) = 1;
                        break;
                    end
                end


                if (foundHard == 0)
                    if (G_lastSid == 0)  % init case.
                        numChoose = sum(chooseY)
                        chooseIdx = randperm(numChoose, 1)
                        count=0;
                        for sid=1:m3_ids
                            if (chooseY(sid) == 1)
                                count++;

                                if (count == chooseIdx)
                                    %specifyId = sid
                                    chooseY = zeros(m3_ids);
                                    chooseY(sid) = 1;
                                    G_lastSid = sid;
                                    break;
                                end
                            end
                        end
                    else
    
%                        rnum = rand(1);
%                        if (rnum>0.66)
%                            chooseY    =  YEstVector == (m3_ids); % default select the champion.
%                        elseif (rnum>0.33)
%                            chooseY    =  YEstVector == (m3_ids - 1);
%                        else
%                            chooseY    =  YEstVector == (m3_ids - 2);
%                        end
    
                        chooseY    =  YEstVector == (m3_ids );
                        [a b] = sort(chooseY);
                        G_lastSid = b(m3_ids);
                    end

                end
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% cxw policy random:
            %%    select more better.
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (policy == 7)

                    chooseY    =  YEstVector == (m3_ids - 0);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%% note
                    %chooseY    =  YEstVector == (m3_ids - 5);
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% cxw policy :
            %%    according by testZone info.
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (policy == 8)
%                global G_testZone;
%                G_testZone(ap, Modelx*3+1) = Modelx;
%                G_testZone(ap, Modelx*3+2) = ROE;
%                        chooseY    =  YEstVector == (m3_ids - 0);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%% note
                    %chooseY    =  YEstVector == (m3_ids - 5);
            end




            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% cxw policy 
            %% pure random.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (policy == -1)
                chooseY=repmat([0], m3_ids, 1);
                randii = 1 + floor(rand(1) * m3_ids);
                chooseY(randii) = 1;
            end




            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% post selection.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            chooseSid = -1;
            minError = 1000;
            GlobalErr = ones(1,m3_ids)*10000;
            for sid=1:m3_ids
                if (chooseY(sid) == 0)
                    continue;
                end

                %% don't consider the invalid one.
                if (YEstVector(sid) <= 1) 
                    continue;
                end

                HError = 0;
                for iPrePreAp=1:2
                    PrePreAp = ap - 1 - iPrePreAp;   %% important!!!
                    indexYest =  m3_header_width + m3_id_width*(sid-1) + m3_index_strategy2_y_estimate;
                    indexYRel =  indexYest - 1;
                    HError(iPrePreAp) = abs(m3(PrePreAp, indexYest) - m3(PrePreAp, indexYRel));
                end

                delta = HError(1) - HError(2);
                GlobalErr(sid) = delta;
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% save to file.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (1)
                [a i] = sort(GlobalErr);  %% Only one item is small than 10000.
                chooseSid = i(specifyRank);
                chooseY(1:m3_ids)  = 0;
                chooseY(chooseSid) = 1;
                
                % *****************************************
                % the ROI module will use this info.
                % in other word, it is output.
                % *****************************************
                m3(ap, m3_index_r1) = chooseSid;

                
                %% formal case, will wirte sid to file.
                %if (policy == 4)
                   % m3(ap, m3_index_r1) =  m3(ap, m3_index_r1) + stockIdStr(chooseSid) / 1000000;
                %end
            end
        end
            
    end  % endof iUp
    end  % endof iLow
end % endof specifyRank


%stockSelect

endfunction

######################################################## FILE.
function  [] = M3_testClear()
######################################################## FILE.
global m3 ;
global m3_ids ;
global nGroups;
global m2_groupSize;

global m3_header_width;
global m3_index_strategy2_y_estimate;
global m3_id_width;

global m3_index_date      ;% = 2;    % header area
global m3_index_envStress ;%= 3;    % header area
global m3_index_r1        ;%= 5;   % header area
global m3_index_r2        ;%= 7;   % header area
global m3_index_r2_roi    ;%= 8;   % header area
global m3_index_r3        ;%= 9;   % header area
global m3_index_r3_roi    ;%= 10;  % header area




global m2;
global m2_header_width;
global m2_index_p;
global m2_id_width;


global G_ROI ;
global stockIdStr ;


m3(1:nGroups, m3_index_r1) = 0;

endfunction





######################################################## FILE.
function  [] = M3_makeRoiTest()
######################################################## FILE.
global m3 ;
global m3_ids ;
global nGroups;
global m2_groupSize;

global m3_header_width;
global m3_index_strategy2_y_estimate;
global m3_id_width;

global m3_index_date      ;% = 2;    % header area
global m3_index_envStress ;%= 3;    % header area
global m3_index_r1        ;%= 5;   % header area

global m3_index_strategy1_roi ;


global m2;
global m2_header_width;
global m2_index_p;
global m2_id_width;


global G_ROI ;
global stockIdStr ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% 
%%  lastAp:    Y_sid... Yest... 
%% 
%%  ap:        Praw
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TESTINFO = {};
testCount = 0;

diffRoi = 0;

for specifyRank=1:1                % super parameter

    ROI_table = 0;
    AP_SID = 0;
    for iLow =1:1                      % super parameter, always be same as std of neighbor.
    for iUp=1:1                        % super parameter
        optLow = -1 * iLow * 0.0025;
        %optUp  = iUp * 0.45 * optLow;
        %optUp  =  1 * iUp  * 0.0045 ;

        validApCount = 0;
        marketInitial = 0;
        marketEnddddd = 0;
        roi = 0;
        TESTCN = 1;

        for testId=1:TESTCN % repeat test.
        for ap=1:nGroups    % very important to howto understand!
                            % ap is the peroid for testify market data.


            thisAPidx  = (m2(:,2) == ap);
            for _i=1:length(thisAPidx)
                if (thisAPidx(_i) == 1)
                    break;
                end
            end
            thisAPidxFirst = _i;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%  Validation!
            %%    reference to chooseSid function.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lastAp = ap - m2_groupSize + 1;           % Note: key point!
            lastAp = ap - 1;                          % Note: key point!
            lastAp = ap    ;                          % Note: key point!
            if (lastAp<1)
                continue;
            end
            chooseSid  = m3(lastAp, m3_index_r1);     % Note!!  Just get from m3.
            if (chooseSid == 0)
                continue;
            end



            indexStart =  m3_header_width + m3_index_strategy2_y_estimate; % est
            indexEnd   =  m3_header_width + m3_ids*m3_id_width;
            YEstIdx    =  indexStart : m3_id_width : indexEnd; 
%            YEstVector  =  m3(lastAp, YEstIdx);
            YRealVector =  m3(lastAp, YEstIdx-1); % no sense!

%
%            % *********************************************
%            % use Y-est-mode-x substitude Y-est!
%            %        dynamic select model.
%            %       
%            % *********************************************
%            global G_ModelSelect;
%            global m3_index_strategy2_b0  ;
%            if (G_ModelSelect == 0)
%                indexStart =  m3_header_width + m3_index_strategy2_y_estimate;
%            else
%                indexStart =  m3_header_width + m3_index_strategy2_b0 + G_ModelSelect; % est y2 (model 1)
%            end
%            indexEnd   =  m3_header_width + m3_ids*m3_id_width;
%            YEstIdx_y1     =  indexStart : m3_id_width : indexEnd; 
%            YEstVector_y1  =  m3(lastAp, YEstIdx_y1);
%
%            [a ix] = sort(YEstVector_y1);
%            tmp(ix) = 1:m3_ids;
%            YEstVector_y1  =  tmp;
%            YEstVector     =  YEstVector_y1;  % Note!  At last, get this.




            %%%%%%%%%%%%%%%%%%%%
            %% make market idx.
            %%%%%%%%%%%%%%%%%%%%
            if (marketInitial == 0)
                    startidx  =  m2_header_width + m2_index_p - 2;
                    endidx    =  m2_header_width + m3_ids*m2_id_width;
                    idxPraw   =  m2(thisAPidxFirst, startidx:m2_id_width:endidx);
                    
                    marketInitial = sum(idxPraw);
            end
            if (ap == nGroups)
                    startidx  =  m2_header_width + m2_index_p - 2;
                    endidx    =  m2_header_width + m3_ids*m2_id_width;
                    idxPraw   =  m2(thisAPidxFirst, startidx:m2_id_width:endidx);
                    
                    marketEnddddd = sum(idxPraw);
            end



            %%%%%%%%%%%%%%%%%%%%
            %% get roi.
            %%%%%%%%%%%%%%%%%%%%
            sid = chooseSid;
            if (1)
                %AP_SID(ap, 1) = sid;
                AP_SID(ap, 1) = stockIdStr(sid);

                AP_SID(ap, 3) = YRealVector(sid);
                %AP_SID(ap, 4) = YEstVector(sid);
                %AP_SID(ap, 5) = abs(YRealVector(sid) - YEstVector(sid));
           

                rawPIdx    =  m2_header_width + m2_index_p - 1 + (sid-1)*m2_id_width;
                rawP       =  m2(thisAPidx, rawPIdx);
                samples = length(rawP);

                if (1)
                    indexRoi = m3_header_width + m3_index_strategy1_roi + (sid-1)*m3_id_width;
                    thisRoi  = m3(ap,indexRoi)   ;
                    roi = roi + thisRoi;
                    vv = thisRoi;
                end


                AP_SID(ap, 6) = vv;
                AP_SID(ap, 7) = rawP(floor(samples*0.6));
                diff = vv - rawP(floor(samples*0.6));
                diffRoi = diffRoi + diff;
                %AP_SID(ap, 8) = vv - rawP(floor(samples/2));
                AP_SID(ap, 8) = diffRoi; 
            end

            validApCount ++;
        end
        end % end of testId
        roi = roi/1;
            
        testCount = testCount + 1;
        ROI_table(iUp, iLow) = roi;
        G_ROI = G_ROI + roi;


        avErr = mean( AP_SID(:, 5) );
        TESTINFO{testCount} = {testCount, ROI_table, specifyRank, optLow, avErr,  AP_SID};
    end  % endof iUp
    end  % endof iLow

end % endof specifyRank



%idxPraw 
marketInitial;
marketEnddddd;
marketROI = (marketEnddddd - marketInitial) / marketInitial;

TESTINFO{testCount+1} = {marketInitial,  marketEnddddd,  marketROI };
format long g
TESTINFO(:)
format short

%validResultSid = m3(:,2) > 0
%stockIdStr(m3(validResultSid,2))
endfunction





######################################################### FILE.
# Global area.
######################################################## FILE.

global m2 = 0;
%if ( exist('runOnLine') )
if (1)
    m2 = dlmread("data_m2.edat", ",");
else
    m2_p1 = dlmread("data_m2.edat.part1", ",");
    m2 = m2_p1;        clear m2_p1;
    
    fl = fopen("data_m2.edat.part2", "r");
    va = fread(fl,100);
    fclose(fl);
    if (length(va) > 0)
        m2_p2 = dlmread("data_m2.edat.part2", ",");
        m2 = [m2 ; m2_p2]; clear m2_p2;
    end
end
global m2_groupSize = 2;         % how many days per group.
global m3_bigGroupSize = 2;
global m3_ids = 0;
global bigGr = 0.;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% m2aux.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global m2aux;
global stockBriefIDStr;
global stockBriefIPO;
global stockBriefGBLT;
global stockBriefNetProfit;
m2aux = dlmread("data_m2aux.edat", ",");
stockBriefIDStr = m2aux(:,1);
stockBriefIPO   = m2aux(:,4);
stockBriefGBLT  = m2aux(:,3);
stockBriefNetProfit  = m2aux(:,6);



global m3 ;
global m3_header_width = 12;     % header area
global m3_id_width = 57;         % header area

global m3_index_date1     = 2;    % header area
global m3_index_date2     = 3;    % header area
global m3_index_envStress = 4;    % header area
global m3_index_global_valid  = 5;   % 
global m3_index_r1            = 6;   % 
global m3_index_market          = 7;   % market price. 
global m3_index_market_1        = 8;   % market roi index
global m3_index_market_2        = 9;   % market roi index (avg, use for fobidden flag).
global m3_index_market_3        = 10;  % header area positive rate.



global m3_index_strategy1_indicator1 = 1;
global m3_index_strategy1_indicator2 = 3;
global m3_index_strategy1_indicator3 = 5;
global m3_index_strategy1_indicator4 = 7;
global m3_index_strategy1_indicator1_pvalue = 2;
global m3_index_strategy1_indicator2_pvalue = 4;
global m3_index_strategy1_indicator3_pvalue = 6;
global m3_index_strategy1_indicator4_pvalue = 8;
global m3_index_strategy1_y = 10;
global m3_index_strategy1_y_estimate = 11;
global m3_index_strategy1_y_estimate2 = 12;
global m3_index_strategy1_b0 = 14;
global m3_index_strategy1_b1 = 16;
global m3_index_strategy1_b2 = 18;
global m3_index_strategy1_b0_pvalue = 15;
global m3_index_strategy1_b1_pvalue = 17;
global m3_index_strategy1_b2_pvalue = 19;
global m3_index_strategy1_roi = 20;
global m3_index_sid     = 21; % sid 

global m3_index_strategy2_indicator1 = 23;
global m3_index_strategy2_indicator2 = 25;
global m3_index_strategy2_indicator3 = 27;
global m3_index_strategy2_indicator4 = 29;
global m3_index_strategy2_indicator1_pvalue = 24;
global m3_index_strategy2_indicator2_pvalue = 26;
global m3_index_strategy2_indicator3_pvalue = 28;
global m3_index_strategy2_indicator4_pvalue = 30;
global m3_index_strategy2_y = 31;
global m3_index_strategy2_y_estimate = 32;
global m3_index_strategy2_y_estimate2 = 33; %% use as sigma!
global m3_index_strategy2_b0 = 34;
global m3_index_strategy2_b1 = 35 ;
global m3_index_strategy2_b2 = 36 ;
global m3_index_strategy2_b3 = 37 ;
global m3_index_strategy2_b4 = 38 ;
global m3_index_strategy2_b5 = 39 ;
global m3_index_strategy2_b6 = 40 ;
global m3_index_strategy2_b7 = 42 ;
global m3_index_strategy2_b8 = 43 ;
global m3_index_strategy2_b0_pvalue = 45 ;
global m3_index_strategy2_b1_pvalue = 46 ;
global m3_index_strategy2_b2_pvalue = 47 ;
global m3_index_strategy2_b3_pvalue = 48 ;
global m3_index_strategy2_b4_pvalue = 49 ;
global m3_index_strategy2_b5_pvalue = 50 ;
global m3_index_strategy2_b6_pvalue = 51 ;
global m3_index_strategy2_b7_pvalue = 52 ;
global m3_index_strategy2_b8_pvalue = 53 ;
global m3_index_strategy2_y_perf1 = 54;
global m3_index_strategy2_y_perf2 = 55;
global m3_index_strategy2_y_perf3 = 56;
global m3_index_indicator_valid   = 57;


# sync with m2Info structure.
global m2_header_width  = 5;    % header area
global m2_header_p_index =4;    % header area
global m2_id_width      = 12;   % header area
global m2_index_Ids     = 2 +1; % header area
global m2_index_sid     = 0 +1; % id area
global m2_index_valid   = 1 +1; % id area
global m2_index_p       = 4 +1; % id area
global m2_index_v       = 9 +1; % id area
global m2_index_jbm_tv  = 10+1; % id area

global G_MaxModels = 8;
global A_table;


global nGroups = max(m2(:,2));
thisGrIndex = (m2(:,2) ==1);
m3_ids = m2(thisGrIndex, m2_index_Ids);
m3_ids = m3_ids(1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read from file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m3 = dlmread("data_m3.edat", "\t");

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% code/script area.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nGroups = nGroups - 1;  % temp!
%m3_ids = 10; %%% tmp!!!

global G_ModelSelect;
global G_ModelSelect_toFile;
global G_ROI;
global G_lastSid;
global G_pg;

G_pg = zeros(1, 8);
G_ROI_TABLE = zeros(1,8);
G_ModelSelect_toFile = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  make sid container.
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M3_compositeSids();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make post work for model3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M3_postModel3();



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  First Run.
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M3_postRank();



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% to follow reality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M3_makeR1R2R3();




if (1)
for iMode = [1:G_MaxModels]

    G_ModelSelect = iMode;
    G_ROI_1 = 0;
    G_ROI_2 = 0;
    G_ROI_3 = 0;
    G_ROI_4 = 0;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  ROI TEST. 1
    %% 
    %%   random (rank 1 or 2) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (0)
        G_ROI = 0;
        TESTCNT = 150;
        for iTest=1:TESTCNT
            M3_makeChooseSid(3);
            M3_makeRoiTest();
            M3_testClear();
        end
        G_ROI_1 = G_ROI / TESTCNT
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  ROI TEST. 2
    %%    (with buffer/pool)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (0)
        G_ROI = 0;
        TESTCNT = 10;
        for iTest=1:TESTCNT
            G_lastSid = 0;
            M3_makeChooseSid(6);
            M3_makeRoiTest();
            M3_testClear();
        end
        G_ROI_2 = G_ROI / TESTCNT
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  ROI TEST. 3
    %%    select champion.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (1)
        ROI_TEST = 3
        G_ROI = 0;
        TESTCNT = 1;             % always be 1.
        for iTest=1:TESTCNT
            M3_makeChooseSid(7); % select champion.
            M3_makeRoiTest();
            M3_testClear();
        end
        G_ROI_3 = G_ROI / TESTCNT
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  ROI TEST. 4
    %%    select  a random one.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (0)
    %if (G_ModelSelect == 1)
        G_ROI = 0;
        TESTCNT = 150;
        for iTest=1:TESTCNT
            M3_makeChooseSid(-1); % select  random
            M3_makeRoiTest();
            M3_testClear();
        end
        G_ROI_4 = G_ROI / TESTCNT
    end
    
    
    
    whichCol = G_ModelSelect;
    if (whichCol == 0)
        whichCol = 8;
    end
    G_ROI_TABLE(1, whichCol) = G_ROI_1;
    G_ROI_TABLE(2, whichCol) = G_ROI_2;
    G_ROI_TABLE(3, whichCol) = G_ROI_3;
    G_ROI_TABLE(4, whichCol) = G_ROI_4;
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% statistic display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_table
B_table = [G_pg; 1.11111*ones(size(G_ROI_TABLE)); G_ROI_TABLE]
dlmwrite("B_Table.edat", B_table, "\t", 'precision', '%2.7f');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% stream to file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M3_writeBack();


return;




######################################################## END OF FILE.
######################################################## END OF FILE.
