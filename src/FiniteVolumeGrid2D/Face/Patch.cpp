#include "Patch.h"

Patch::Patch(const std::string &name, Label id, const std::shared_ptr<PatchRegistry> &registry)
    :
      Group(name)
{
    id_ = id;

    if(registry)
        registry_ = registry;
    else
        registry_ = std::shared_ptr<PatchRegistry>(new PatchRegistry());
}

void Patch::add(const Face &face)
{
    auto insertion = registry_->insert(std::make_pair(face.id(), std::ref(*this)));

    if(insertion.second)
        FaceGroup::add(face);
    else
    {
        insertion.first->second.get().remove(face);
        add(face);
    }
}

void Patch::remove(const Face &face)
{
    if(isInGroup(face))
    {
        registry_->erase(face.id());
        FaceGroup::remove(face);
    }
}

void Patch::clear()
{
    for(const Face& face: *this)
        registry_->erase(face.id());

    FaceGroup::clear();
}

void Patch::setRegistry(std::shared_ptr<Patch::PatchRegistry> &registry)
{
    std::vector<Ref<const Face>> items = this->items();
    clear();
    registry_ = registry;
    FaceGroup::add(items.begin(), items.end());
}
