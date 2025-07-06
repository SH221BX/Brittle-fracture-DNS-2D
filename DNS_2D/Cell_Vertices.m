function [all_cell_verts, int_points] = Cell_Vertices(C, V, constraints, bx)
    all_cell_verts = [];
    int_points = [];

    for cell_no = 1:size(C, 1)
        cell_vertex_indices = C{cell_no};
        cell_verts = V(cell_vertex_indices,:);
        cell_verts = cell_verts(~isinf(cell_verts(:,1)), :);
        cell_verts = cell_verts(~isinf(cell_verts(:,2)), :);
        if any(cell_verts(:) < 0) || any(cell_verts(:) > bx)

            [xint, yint] = polyxpoly(cell_verts(:,1), cell_verts(:,2), ...
                constraints(:,1), constraints(:,2));

            int_points = [int_points; xint, yint];
            int_points = unique(int_points, 'rows');
            cell_verts = [cell_verts; xint, yint; int_points];
            cell_verts = unique(cell_verts, 'rows');
        end
        cell_verts = cell_verts(all(cell_verts >= 0 & cell_verts <= bx, 2), :);
        all_cell_verts = [all_cell_verts; cell_verts];
        all_cell_verts = unique(all_cell_verts, 'rows');
    end
end
