#include "Numbered.h"

Numbered::Numbered(Numbered const & that)
{
    this->id = that.id;
}

Numbered::Numbered()
{
    static int id_static = 0;
    this->id = id_static;
    id_static += 1;
}

int Numbered::get_id() const
{
    return id;
}

bool Numbered::operator==(Numbered const & that)
{
    return id == that.id;
}

bool Numbered::operator!=(Numbered const & that)
{
    return !(*this == that);
}

Numbered & Numbered::operator=(Numbered const & that)
{
    id = that.id;
    return *this;
}
