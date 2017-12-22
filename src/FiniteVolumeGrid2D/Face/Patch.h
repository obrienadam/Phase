#ifndef PATCH_H
#define PATCH_H

#include "FaceGroup.h"

class Patch: public FaceGroup
{
public:

    typedef std::unordered_map<Label, Ref<Patch>> PatchRegistry;

    Patch(const std::string &name, Label id, const std::shared_ptr<PatchRegistry> &registry);

    ~Patch()
    {
        clear(); //- Make sure that the registry is updated
    }

    void add(const Face &face);

    void add(const FaceGroup& group);

    void remove(const Face &face);

    void remove(const FaceGroup& faces);

    void clear();

    Label id() const
    { return id_; }

    void setRegistry(std::shared_ptr<PatchRegistry>& registry);

    std::shared_ptr<PatchRegistry> registry()
    { return registry_; }

private:

    Label id_;

    std::shared_ptr<PatchRegistry> registry_;
};

#endif
