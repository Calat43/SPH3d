#pragma once

class Cell;

namespace nl_idic
{
    double sound_speed_asterisk(Cell * cell);

    double lambda_asterisk(Cell * cell);

    double mach_number_asterisk(Cell * cell);
    double knudsen_number_asterisk(Cell * cell);
    double reynolds_number_asterisk(Cell * cell);

    double t_stop_asterisk(Cell * cell);
}