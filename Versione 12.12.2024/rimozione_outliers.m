function [x,y] = rimozione_outliers(trajectory)
    % If in the given trajectory there's a point superimposed with another, this function removes that point
    % trajectory    is the input trajectory
    % 
    % (x, y)        is the output trajectory

    x(1)=trajectory(1,1);
    y(1)=trajectory(1,2);
    j = 2;
    for i=2:1:length(trajectory(:,1))
        if abs(trajectory(i,1)-trajectory(i-1,1))>1e-6 && abs(trajectory(i,2)-trajectory(i-1,2))>1e-6
           x(j)=trajectory(i,1);
           y(j)=trajectory(i,2);
           j = j+1;
        end 
    end
end