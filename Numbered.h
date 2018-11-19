#pragma once

class Numbered
{
public:
    Numbered();
    Numbered(Numbered const & that);
    int get_id() const;

    Numbered & operator=(Numbered const & that);

    bool operator==(Numbered const & that);
    bool operator!=(Numbered const & that);
private:
    int id;
};