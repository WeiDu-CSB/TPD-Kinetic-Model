function [value,isterminal,direction] = events_time(t,y,parameter)

etime=toc;

value = etime-60;
isterminal=1;
direction=0;