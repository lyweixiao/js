function[Distance] = compute_distance1(BS, UE)
Distance = sqrt((BS.location_x-UE.location_x)^2+(BS.location_y-UE.location_y)^2);