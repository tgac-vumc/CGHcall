setClass("cghRaw",
        contains    = "eSet",
        prototype   = prototype(new("VersionedBiobase",
        versions    = c(classVersion("eSet"), cghRaw="1.0.0")))
)

setClass("cghSeg",
        contains    = "eSet",
        prototype   = prototype(new("VersionedBiobase",
        versions    = c(classVersion("eSet"), cghSeg="1.0.0")))
)

setClass("cghCall",
        contains    = "eSet",
        prototype   = prototype(new("VersionedBiobase",
        versions    = c(classVersion("eSet"), cghCall="1.0.0")))
)
